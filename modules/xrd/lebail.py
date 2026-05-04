"""
modules/xrd/lebail.py
Le Bail profile fitting engine — pymatgen-powered.

Uses pymatgen's XRDCalculator for correct structure factor + LP computation
to seed I_hkl. Falls back to our crystallography.py if pymatgen unavailable.

Outer loop:
  A) Le Bail I_hkl update  (intensity partitioning)
  B) Global scale factor   (linear regression)
  C) Background update     (Chebyshev polynomial)
  D) Cell + profile params (scipy least_squares)
"""

import math, warnings
import numpy as np
from scipy.optimize import least_squares

from .crystallography import (
    generate_reflections, chebyshev_background,
    caglioti_fwhm, tch_fwhm_eta, scherrer_size, size_from_Y,
    pseudo_voigt,
    compute_fit_statistics, d_spacing, cell_volume,
    generate_reflections_rietveld, compute_rietveld_intensities,
    structure_factor_sq_dw, parse_cif as _parse_cif_cryst,
    molar_mass_from_formula, expand_sites_from_cif,
)

# ─────────────────────────────────────────────────────────────────────────────
# PYMATGEN INTEGRATION
# ─────────────────────────────────────────────────────────────────────────────

def _try_import_pymatgen():
    try:
        from pymatgen.core import Structure
        from pymatgen.analysis.diffraction.xrd import XRDCalculator
        return Structure, XRDCalculator
    except ImportError:
        return None, None


def _cif_to_structure(cif_text):
    """Parse CIF text into a pymatgen Structure. Returns None on failure."""
    Structure, _ = _try_import_pymatgen()
    if Structure is None:
        return None
    import tempfile, os
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif',
                                         delete=False) as f:
            f.write(cif_text)
            path = f.name
        struct = Structure.from_file(path)
        os.unlink(path)
        return struct
    except Exception:
        return None


def get_pymatgen_intensities(cif_text, wavelength, tt_min, tt_max):
    """
    Use pymatgen XRDCalculator to get per-reflection intensities
    with proper structure factors and LP correction.

    Returns dict: {two_theta_rounded: intensity} or None if unavailable.
    """
    _, XRDCalculator = _try_import_pymatgen()
    if XRDCalculator is None:
        return None

    struct = _cif_to_structure(cif_text)
    if struct is None:
        return None

    try:
        calc    = XRDCalculator(wavelength=wavelength)
        pattern = calc.get_pattern(struct, scaled=False,
                                    two_theta_range=(tt_min, tt_max))
        # pattern.x = 2theta positions, pattern.y = intensities
        intensity_map = {}
        for tt, I in zip(pattern.x, pattern.y):
            intensity_map[round(float(tt), 3)] = float(I)
        return intensity_map
    except Exception:
        return None


def seed_I_hkl_from_pymatgen(refs, intensity_map, tt_r, y_r, bg_est):
    """
    Seed I_hkl values using observed peak heights, but informed by pymatgen:
    - Peaks that pymatgen says have essentially zero intensity are zeroed out
    - All other peaks are seeded from the local observed intensity at that 2θ
    
    This gives Le Bail a physically informed starting point without
    over-constraining the intensity ratios (which Le Bail should determine
    from the data, not from theory).
    """
    I_hkl = np.zeros(len(refs))
    pm_max = max(intensity_map.values()) if intensity_map else 1.0

    for k, (tt_peak, d, hkl, mult) in enumerate(refs):
        # Check if pymatgen says this reflection exists
        best_I = None
        best_dist = 0.15
        for pm_tt, pm_I in (intensity_map or {}).items():
            dist = abs(pm_tt - tt_peak)
            if dist < best_dist:
                best_dist = dist
                best_I = pm_I

        if best_I is not None and best_I / pm_max < 0.001:
            # Pymatgen says this peak is essentially zero — keep it near zero
            I_hkl[k] = 0.01
        else:
            # Seed from local observed data height at this peak position
            near = np.abs(tt_r - tt_peak) < 0.5
            if near.any():
                I_hkl[k] = max(y_r[near].max() - bg_est, 1.0)
            else:
                I_hkl[k] = 1.0

    return np.maximum(I_hkl, 0.01)


# ─────────────────────────────────────────────────────────────────────────────
# PROFILE FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# PROFILE FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def _get_profiles(tt_arr, refs, U, V, W, eta, zero=0.0, X=0.0, Y=0.0,
                   window_factor=15.0):
    """
    Unit-normalised pseudo-Voigt profiles, one per reflection.

    If X or Y are non-zero, uses Thompson-Cox-Hastings (TCH) model where
    eta is computed per-peak from Gaussian (U,V,W) and Lorentzian (X,Y)
    widths.  Otherwise uses fixed-eta Caglioti pseudo-Voigt.

    window_factor : float
        Profile evaluation window as a multiple of FWHM.  Use 15.0 (default)
        for refinement (numerical accuracy) and 3.0 for display patterns
        (suppresses Lorentzian tail cross-contamination between phases).
    """
    if not refs:
        return []

    use_tch = (X != 0.0 or Y != 0.0)
    tt_s  = tt_arr - zero
    n_pts = len(tt_arr)
    n_ref = len(refs)
    profiles = [None] * n_ref

    if not use_tch:
        eta = float(np.clip(eta, 0, 1))

    for k, (tt_p, d, hkl, mult) in enumerate(refs):
        if use_tch:
            fwhm, eta_k = tch_fwhm_eta(tt_p, U, V, W, X, Y)
        else:
            fwhm  = max(caglioti_fwhm(tt_p, U, V, W), 0.005)
            eta_k = eta

        sigma_g = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        gamma_l = fwhm / 2.0
        window  = window_factor * fwhm

        prof = np.zeros(n_pts)
        msk  = np.abs(tt_s - tt_p) < window
        if msk.any():
            dx      = tt_s[msk] - tt_p
            gauss   = np.exp(-0.5 * (dx / sigma_g) ** 2)
            lor     = 1.0 / (1.0 + (dx / gamma_l) ** 2)
            prof[msk] = eta_k * lor + (1.0 - eta_k) * gauss
        profiles[k] = prof

    return profiles


# ─────────────────────────────────────────────────────────────────────────────
# CELL PARAMETER HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _cell_free(phase):
    """Return (free_values, free_names) for this crystal system."""
    a  = float(phase['a'])
    b  = float(phase.get('b') or a)
    c  = float(phase.get('c') or a)
    be = float(phase.get('beta', 90.0)  or 90.0)
    al = float(phase.get('alpha', 90.0) or 90.0)
    ga = float(phase.get('gamma', 90.0) or 90.0)
    sys_ = (phase.get('system') or 'triclinic').lower()
    if sys_ == 'cubic':
        return [a], ['a']
    elif sys_ == 'tetragonal':
        return [a, c], ['a', 'c']
    elif sys_ in ('hexagonal', 'trigonal'):
        return [a, c], ['a', 'c']
    elif sys_ == 'orthorhombic':
        return [a, b, c], ['a', 'b', 'c']
    elif sys_ == 'monoclinic':
        return [a, b, c, be], ['a', 'b', 'c', 'beta']
    else:
        return [a, b, c, al, be, ga], ['a', 'b', 'c', 'alpha', 'beta', 'gamma']


def _full_cell(free_vals, free_names, phase):
    """Reconstruct (a,b,c,al,be,ga) from free params + symmetry constraints."""
    d = {
        'a': float(phase['a']), 'b': float(phase.get('b') or phase['a']),
        'c': float(phase.get('c') or phase['a']),
        'alpha': float(phase.get('alpha', 90.0) or 90.0),
        'beta':  float(phase.get('beta',  90.0) or 90.0),
        'gamma': float(phase.get('gamma', 90.0) or 90.0),
    }
    for name, val in zip(free_names, free_vals):
        d[name] = val
    sys_ = (phase.get('system') or 'triclinic').lower()
    # Apply all symmetry constraints — angles as well as lengths
    if sys_ == 'cubic':
        d['b'] = d['a']; d['c'] = d['a']
        d['alpha'] = d['beta'] = d['gamma'] = 90.0
    elif sys_ == 'tetragonal':
        d['b'] = d['a']
        d['alpha'] = d['beta'] = d['gamma'] = 90.0
    elif sys_ == 'hexagonal':
        d['b'] = d['a']
        d['alpha'] = d['beta'] = 90.0; d['gamma'] = 120.0
    elif sys_ == 'trigonal':
        d['b'] = d['a']
        d['alpha'] = d['beta'] = 90.0; d['gamma'] = 120.0
    elif sys_ == 'orthorhombic':
        d['alpha'] = d['beta'] = d['gamma'] = 90.0
    elif sys_ == 'monoclinic':
        d['alpha'] = d['gamma'] = 90.0  # only beta is free
    return d['a'], d['b'], d['c'], d['alpha'], d['beta'], d['gamma']


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _filter_tick_positions(refs, I_hkl, threshold_frac=1e-3):
    """Filter tick positions, keeping only reflections with significant intensity.

    Removes ghost peaks whose refined I_hkl is < 0.1% of the strongest peak.
    Works for both Le Bail refs (list of tuples) and Rietveld refs (list of dicts).
    """
    if len(refs) == 0:
        return []
    max_I = max(I_hkl) if len(I_hkl) > 0 else 0
    threshold = max_I * threshold_frac if max_I > 0 else 0
    ticks = []
    for k, ref in enumerate(refs):
        if k < len(I_hkl) and I_hkl[k] > threshold:
            tt_val = ref['two_theta'] if isinstance(ref, dict) else ref[0]
            ticks.append(round(tt_val, 3))
    return ticks


# ─────────────────────────────────────────────────────────────────────────────
# MAIN REFINEMENT
# ─────────────────────────────────────────────────────────────────────────────

def run_lebail(tt, y_obs, sigma, phases, wavelength,
               tt_min=None, tt_max=None,
               n_bg_coeffs=6, max_outer=15,
               progress_callback=None):
    """
    Run Le Bail refinement.

    Parameters
    ----------
    tt, y_obs, sigma : np.ndarray — data arrays
    phases   : list of dicts with keys:
                 name, system, spacegroup_number, a, [b, c, alpha, beta, gamma]
                 optional: cif_text (enables pymatgen seeding)
    wavelength : float — Å
    tt_min/max : float — refinement range
    n_bg_coeffs: int
    max_outer  : int

    Returns dict with tt, y_obs, y_calc, y_background, phase_patterns,
                      residuals, statistics, phase_results, zero_shift
    """
    if tt_min is None: tt_min = float(tt.min())
    if tt_max is None: tt_max = float(tt.max())

    mask   = (tt >= tt_min) & (tt <= tt_max)
    tt_r   = tt[mask]
    y_r    = y_obs[mask]
    sig_r  = (sigma[mask] if sigma is not None
              else np.sqrt(np.maximum(y_r, 1.0)))
    weights = 1.0 / np.maximum(sig_r**2, 1e-6)

    if progress_callback: progress_callback('Generating reflections...')

    # ── Per-phase state ───────────────────────────────────────────────────
    phase_state = []
    bg_est = float(np.percentile(y_r, 5))

    for ph in phases:
        free_v, free_n = _cell_free(ph)
        a, b, c, al, be, ga = _full_cell(free_v, free_n, ph)

        # Extract atom sites from CIF text if available
        # (enables structure-factor-based reflection filtering).
        # Use pymatgen expansion to get full unit cell for correct F²,
        # then fall back to raw parse_cif (asymmetric unit only).
        sites = ph.get('sites')
        if not sites and ph.get('cif_text'):
            sites = expand_sites_from_cif(ph['cif_text'])
            if not sites:
                try:
                    _parsed = _parse_cif_cryst(ph['cif_text'])
                    sites = _parsed.get('sites') or None
                except Exception:
                    sites = None
        ph['_sites'] = sites or None  # stash for cell-refinement regeneration

        refs = generate_reflections(
            a, b, c, al, be, ga,
            ph.get('system', 'triclinic'),
            ph.get('spacegroup_number', 1),
            wavelength, tt_min, tt_max, hkl_max=12,
            sites=ph['_sites']
        )

        # Seed I_hkl — use pymatgen if CIF is available
        cif_text = ph.get('cif_text')
        pm_map   = None
        if cif_text:
            if progress_callback:
                progress_callback(f"Computing structure factors for {ph.get('name','?')}...")
            pm_map = get_pymatgen_intensities(cif_text, wavelength, tt_min, tt_max)

        if pm_map:
            I_hkl = seed_I_hkl_from_pymatgen(refs, pm_map, tt_r, y_r, bg_est)
            seeded = 'pymatgen'
        else:
            # Fallback: use observed peak heights
            I_hkl = np.zeros(len(refs))
            for k, (tt_p, d, hkl, mult) in enumerate(refs):
                near = np.abs(tt_r - tt_p) < 0.5
                I_hkl[k] = max(y_r[near].max() - bg_est, 1.0) if near.any() else 1.0
            seeded = 'observed'

        phase_state.append({
            'ph':     ph,
            'refs':   refs,
            'I_hkl':  I_hkl,
            'S':  1.0,  # per-phase scale factor
            'U':  ph.get('U_init',   0.01),
            'V':  ph.get('V_init',  -0.01),
            'W':  ph.get('W_init',   0.15),
            'X':  ph.get('X_init',   0.0),   # Lorentzian strain
            'Y':  ph.get('Y_init',   0.1),   # Lorentzian size (Scherrer)
            'eta':ph.get('eta_init', 0.5),
            'free_v': free_v,
            'free_n': free_n,
            'seeded': seeded,
        })

    n_bg  = n_bg_coeffs
    bg_c  = np.zeros(n_bg); bg_c[0] = bg_est
    zero  = 0.0

    def total_calc(states, bg_c_, zero_):
        bg_  = np.maximum(chebyshev_background(tt_r, bg_c_, tt_min, tt_max), 0)
        pat_ = np.zeros(len(tt_r))
        for st in states:
            profs = _get_profiles(tt_r, st['refs'],
                                   st['U'], st['V'], st['W'],
                                   st.get('eta', 0.5),
                                   zero_, st['X'], st['Y'])
            for k in range(len(st['refs'])):
                pat_ += st['S'] * st['I_hkl'][k] * profs[k]
        return pat_ + bg_, bg_, pat_

    def rwp(yo, yc, w):
        d = yo - yc
        return math.sqrt(np.sum(w*d**2) / np.sum(w*yo**2)) * 100

    if progress_callback: progress_callback('Running Le Bail iterations...')

    prev_rwp = 999.0
    for outer in range(max_outer):

        # ── A: Le Bail I_hkl update ───────────────────────────────────────
        # Cache profiles once per outer iteration (major speedup)
        all_profs = [_get_profiles(tt_r, st['refs'],
                                    st['U'], st['V'], st['W'],
                                    st.get('eta', 0.5), zero,
                                    st['X'], st['Y'])
                     for st in phase_state]

        for inner in range(60):
            bg_cur  = np.maximum(chebyshev_background(tt_r, bg_c, tt_min, tt_max), 0)
            pat_all = np.zeros(len(tt_r))
            for i_ph, st in enumerate(phase_state):
                for k in range(len(st['refs'])):
                    pat_all += st['S'] * st['I_hkl'][k] * all_profs[i_ph][k]
            sum_c  = np.maximum(pat_all, 1e-8)
            y_nobg = np.maximum(y_r - bg_cur, 0.1)

            max_rel_chg = 0.0
            for i_ph, st in enumerate(phase_state):
                I_new = np.zeros(len(st['refs']))
                for k in range(len(st['refs'])):
                    phi_k = all_profs[i_ph][k]
                    if phi_k.sum() < 1e-10:
                        I_new[k] = st['I_hkl'][k]; continue
                    ratio  = st['S'] * st['I_hkl'][k] * phi_k / sum_c
                    numer  = np.sum(weights * y_nobg * ratio)
                    denom  = np.sum(weights * phi_k) * st['S']
                    I_new[k] = max(numer / max(denom, 1e-10), 1e-6)
                I_new = np.clip(I_new, 1e-6, 1e7)
                # Relative convergence check (scale-independent)
                rel_chg = np.max(np.abs(I_new - st['I_hkl']) /
                                 np.maximum(st['I_hkl'], 1e-6))
                max_rel_chg = max(max_rel_chg, rel_chg)
                st['I_hkl'] = I_new
            if max_rel_chg < 0.005:  # 0.5% relative change threshold
                break

        # ── B: Update per-phase scale factors ─────────────────────────────
        bg_cur = np.maximum(chebyshev_background(tt_r, bg_c, tt_min, tt_max), 0)
        y_nobg2 = y_r - bg_cur
        for i_ph, st in enumerate(phase_state):
            profs = _get_profiles(tt_r, st['refs'],
                                   st['U'], st['V'], st['W'],
                                   st.get('eta', 0.5), zero,
                                   st['X'], st['Y'])
            # Other phases' contribution (fixed for this phase's scale update)
            other_pat = np.zeros(len(tt_r))
            for j, st2 in enumerate(phase_state):
                if j == i_ph:
                    continue
                p2 = _get_profiles(tt_r, st2['refs'],
                                    st2['U'], st2['V'], st2['W'],
                                    st2.get('eta', 0.5), zero,
                                    st2['X'], st2['Y'])
                for k in range(len(st2['refs'])):
                    other_pat += st2['S'] * st2['I_hkl'][k] * p2[k]
            # This phase's unscaled pattern
            this_pat = np.zeros(len(tt_r))
            for k in range(len(st['refs'])):
                this_pat += st['I_hkl'][k] * profs[k]
            # Solve: S_i = argmin Σ w*(y - bg - other - S_i*this)²
            resid_target = y_nobg2 - other_pat
            S_num = np.sum(weights * resid_target * this_pat)
            S_den = np.sum(weights * this_pat**2)
            st['S'] = max(S_num / max(S_den, 1e-10), 1e-6)

        # ── C: Update background ──────────────────────────────────────────
        pat_total = np.zeros(len(tt_r))
        for st in phase_state:
            profs = _get_profiles(tt_r, st['refs'],
                                   st['U'], st['V'], st['W'],
                                   st.get('eta', 0.5), zero,
                                   st['X'], st['Y'])
            for k in range(len(st['refs'])):
                pat_total += st['S'] * st['I_hkl'][k] * profs[k]
        def resid_bg(x):
            bg_ = np.maximum(chebyshev_background(tt_r, x, tt_min, tt_max), 0)
            return (y_r - pat_total - bg_) * np.sqrt(weights)
        rb = least_squares(resid_bg, bg_c, method='trf',
                           max_nfev=200, ftol=1e-10, verbose=0)
        bg_c = rb.x

        # ── D: SIMULTANEOUS cell + profile refinement (all phases) ──────
        if progress_callback and outer == 0:
            progress_callback('Refining cell parameters and peak shape...')

        bg_fixed = np.maximum(chebyshev_background(tt_r, bg_c, tt_min, tt_max), 0)

        # Build joint parameter vector for all phases:
        # Per phase: [cell_params..., S, U, V, W, X, Y]
        # Global:    [zero]
        _ref_caches_lb = [{'fv': None, 'refs': st['refs']} for st in phase_state]
        phase_layouts_lb = []
        x0_lb, lo_lb, hi_lb = [], [], []

        for i_ph, st in enumerate(phase_state):
            ph = st['ph']
            free_n = st['free_n']
            offset = len(x0_lb)

            fv0 = list(st['free_v'])
            n_cell = len(fv0)
            x0_lb.extend(fv0)
            lo_lb.extend([v * 0.94 for v in fv0])
            hi_lb.extend([v * 1.06 for v in fv0])
            for i_p, nm in enumerate(free_n):
                if nm in ('alpha', 'beta', 'gamma'):
                    lo_lb[offset + i_p] = fv0[i_p] - 5
                    hi_lb[offset + i_p] = fv0[i_p] + 5

            # S, U, V, W, X, Y
            x0_lb.extend([st['S'], st['U'], st['V'], st['W'], st['X'], st['Y']])
            lo_lb.extend([1e-6, 0.0, -1.0, 0.005, 0.0, 0.0])
            hi_lb.extend([1e4,  5.0,  1.0, 3.0,   2.0, 5.0])

            phase_layouts_lb.append((offset, n_cell))

        # Global zero
        x0_lb.append(zero); lo_lb.append(-0.3); hi_lb.append(0.3)

        def resid_lb_joint(x):
            zero_ = x[-1]
            pat = np.zeros(len(tt_r))

            for i_ph, st in enumerate(phase_state):
                off, n_cell = phase_layouts_lb[i_ph]
                ph = st['ph']
                free_n = st['free_n']
                sys_ = (ph.get('system') or 'triclinic').lower()

                fv = x[off:off + n_cell]
                idx = off + n_cell
                S_ph = max(x[idx], 1e-6)
                U_ = x[idx+1]; V_ = x[idx+2]; W_ = max(x[idx+3], 0.005)
                X_ = x[idx+4]; Y_ = x[idx+5]

                cache = _ref_caches_lb[i_ph]
                if (cache['fv'] is None or
                        np.max(np.abs(np.array(fv) - np.array(cache['fv']))) > 1e-4):
                    a_,b_,c_,al_,be_,ga_ = _full_cell(fv, free_n, ph)
                    try:
                        cache['refs'] = generate_reflections(
                            a_,b_,c_,al_,be_,ga_, sys_,
                            ph.get('spacegroup_number', 1),
                            wavelength, tt_min, tt_max, hkl_max=12,
                            sites=ph.get('_sites'))
                    except Exception:
                        pass
                    cache['fv'] = list(fv)

                refs_l = cache['refs']
                n = min(len(refs_l), len(st['I_hkl']))
                profs_l = _get_profiles(tt_r, refs_l, U_, V_, W_, 0.5,
                                        zero_, X_, Y_)
                for k in range(n):
                    pat += S_ph * st['I_hkl'][k] * profs_l[k]

            return (y_r - pat - bg_fixed) * np.sqrt(weights)

        try:
            rs = least_squares(resid_lb_joint, x0_lb, bounds=(lo_lb, hi_lb),
                               method='trf', max_nfev=500,
                               ftol=1e-6, xtol=1e-6, verbose=0)
            x_out = rs.x
        except Exception:
            x_out = np.array(x0_lb)

        # Unpack results back into phase_state
        zero = x_out[-1]
        for i_ph, st in enumerate(phase_state):
            off, n_cell = phase_layouts_lb[i_ph]
            ph = st['ph']
            free_n = st['free_n']
            sys_ = (ph.get('system') or 'triclinic').lower()

            fv_new = x_out[off:off + n_cell]
            idx = off + n_cell
            S_n = max(x_out[idx], 1e-6)
            U_n = x_out[idx+1]; V_n = x_out[idx+2]
            W_n = max(x_out[idx+3], 0.005)
            X_n = x_out[idx+4]; Y_n = x_out[idx+5]

            a_n, b_n, c_n, al_n, be_n, ga_n = _full_cell(fv_new, free_n, ph)
            try:
                refs_new = generate_reflections(
                    a_n, b_n, c_n, al_n, be_n, ga_n, sys_,
                    ph.get('spacegroup_number', 1),
                    wavelength, tt_min, tt_max, hkl_max=12,
                    sites=ph.get('_sites'))
            except Exception:
                refs_new = st['refs']

            n_new = len(refs_new)
            I_new = np.ones(n_new)
            n_copy = min(n_new, len(st['I_hkl']))
            I_new[:n_copy] = st['I_hkl'][:n_copy]

            for nm, vl in zip(free_n, fv_new):
                ph[nm] = vl
            if sys_ in ('cubic',):
                ph['b'] = ph['c'] = a_n
            elif sys_ in ('hexagonal', 'trigonal', 'tetragonal'):
                ph['b'] = a_n
            ph['c'] = c_n

            st.update({'free_v': fv_new, 'refs': refs_new, 'I_hkl': I_new,
                       'S': S_n, 'U': U_n, 'V': V_n, 'W': W_n,
                       'X': X_n, 'Y': Y_n})

        # ── E: Quick I_hkl re-partitioning after cell refinement ─────────
        # Peaks may have shifted; re-update intensities for consistency.
        _profs_e = [_get_profiles(tt_r, st['refs'],
                                   st['U'], st['V'], st['W'], st.get('eta', 0.5),
                                   zero, st['X'], st['Y'])
                    for st in phase_state]
        bg_e = np.maximum(chebyshev_background(tt_r, bg_c, tt_min, tt_max), 0)
        for _inner_e in range(10):
            pat_e = np.zeros(len(tt_r))
            for i_ph, st in enumerate(phase_state):
                for k in range(len(st['refs'])):
                    pat_e += st['S'] * st['I_hkl'][k] * _profs_e[i_ph][k]
            sum_e = np.maximum(pat_e, 1e-8)
            y_nob = np.maximum(y_r - bg_e, 0.1)
            max_chg = 0.0
            for i_ph, st in enumerate(phase_state):
                I_up = np.zeros(len(st['refs']))
                for k in range(len(st['refs'])):
                    phi = _profs_e[i_ph][k]
                    if phi.sum() < 1e-10:
                        I_up[k] = st['I_hkl'][k]; continue
                    ratio = st['S'] * st['I_hkl'][k] * phi / sum_e
                    I_up[k] = max(np.sum(weights * y_nob * ratio) /
                                  max(np.sum(weights * phi) * st['S'], 1e-10), 1e-6)
                I_up = np.clip(I_up, 1e-6, 1e7)
                max_chg = max(max_chg, np.max(np.abs(I_up - st['I_hkl']) /
                              np.maximum(st['I_hkl'], 1e-6)))
                st['I_hkl'] = I_up
            if max_chg < 0.005:
                break

        yc_cur, _, _ = total_calc(phase_state, bg_c, zero)
        cur_rwp = rwp(y_r, yc_cur, weights)
        if progress_callback:
            progress_callback(f'Outer iteration {outer+1}: Rwp = {cur_rwp:.2f}%')
        if abs(prev_rwp - cur_rwp) < 0.1:  # stop when Rwp changes less than 0.1%
            break
        prev_rwp = cur_rwp

    # ── Final results ─────────────────────────────────────────────────────
    if progress_callback: progress_callback('Computing statistics...')

    yc_f, bg_f, pat_f_total = total_calc(phase_state, bg_c, zero)
    diff_f = y_r - yc_f
    # n_params: cell params + profile params + bg.  Le Bail I_hkl are
    # not counted — they are constrained by profile overlap and do not
    # behave as independent free parameters (see David 2004).
    n_params = (sum(len(st['free_n']) + 6  # +6: S, U, V, W, X, Y
                    for st in phase_state) + n_bg + 1)  # +1: zero
    stats = compute_fit_statistics(y_r, yc_f, weights, n_params)

    phase_patterns = []
    phase_results  = []

    # ── Weight fractions via Hill & Howard (1987) ────────────────────────
    # W_α = S_α · Z_α · M_α · V_α  /  Σ_i(S_i · Z_i · M_i · V_i)
    # If Z or M are unavailable for any phase, fall back to S·V² (area-based).
    zmv_values = []
    use_zmv = True
    for st in phase_state:
        ph = st['ph']
        a_ = float(ph.get('a', 4.0))
        b_ = float(ph.get('b') or a_)
        c_ = float(ph.get('c') or a_)
        al_ = float(ph.get('alpha', 90.0) or 90.0)
        be_ = float(ph.get('beta',  90.0) or 90.0)
        ga_ = float(ph.get('gamma', 90.0) or 90.0)
        V = cell_volume(a_, b_, c_, al_, be_, ga_)
        Z = ph.get('Z')
        M = molar_mass_from_formula(ph.get('formula', ''))
        if Z and M:
            zmv_values.append(float(st['S']) * float(Z) * float(M) * V)
        else:
            use_zmv = False
            break

    if not use_zmv:
        # Fallback: use integrated pattern areas (semi-quantitative)
        zmv_values = []
        for st in phase_state:
            profs_tmp = _get_profiles(tt_r, st['refs'],
                                       st['U'], st['V'], st['W'],
                                       st.get('eta', 0.5), zero,
                                       st['X'], st['Y'])
            area = st['S'] * sum(st['I_hkl'][k] * profs_tmp[k].sum()
                                 for k in range(len(st['refs'])))
            zmv_values.append(float(area))

    total_zmv = sum(zmv_values) or 1e-10

    for i_ph, st in enumerate(phase_state):
        # Use tight profiles (3× FWHM) for display to prevent Lorentzian
        # tails from one phase bleeding into another phase's peak regions.
        display_profs = _get_profiles(tt_r, st['refs'],
                                       st['U'], st['V'], st['W'],
                                       st.get('eta', 0.5), zero,
                                       st['X'], st['Y'],
                                       window_factor=3.0)
        pat_ph = st['S'] * sum(st['I_hkl'][k]*display_profs[k]
                               for k in range(len(st['refs'])))
        phase_patterns.append(pat_ph.tolist())

        ph = st['ph']
        a  = float(ph.get('a', 4.0))
        b  = float(ph.get('b') or a)
        c  = float(ph.get('c') or a)
        al = float(ph.get('alpha', 90.0) or 90.0)
        be = float(ph.get('beta',  90.0) or 90.0)
        ga = float(ph.get('gamma', 90.0) or 90.0)

        strongest_tt = (max(st['refs'], key=lambda r: r[3])[0]
                        if st['refs'] else 39.0)
        # Crystallite size: prefer TCH Y parameter (direct size extraction),
        # fall back to total Caglioti FWHM if Y is negligible.
        Y_val = st['Y']
        if Y_val > 0.01:
            cryst_A = size_from_Y(Y_val, wavelength)
        else:
            fwhm_tot, _ = tch_fwhm_eta(strongest_tt, st['U'], st['V'],
                                         st['W'], st['X'], st['Y'])
            cryst_A = scherrer_size(fwhm_tot, strongest_tt, wavelength)

        fwhm_sc, eta_sc = tch_fwhm_eta(strongest_tt, st['U'], st['V'],
                                         st['W'], st['X'], st['Y'])

        wt_frac = (zmv_values[i_ph] / total_zmv) * 100

        _sg_sym_r  = (ph.get('spacegroup', '') or
                      f"SG{ph.get('spacegroup_number', 1)}")
        _base_r    = ph.get('name') or ph.get('formula') or f'Phase {i_ph+1}'
        if _sg_sym_r and _sg_sym_r.replace(' ','') not in _base_r.replace(' ',''):
            _display_r = f"{_base_r} {_sg_sym_r}"
        else:
            _display_r = _base_r

        phase_results.append({
            'name':              _display_r,
            'cod_id':            ph.get('cod_id', ''),
            'formula':           ph.get('formula', ''),
            'a': round(a,5), 'b': round(b,5), 'c': round(c,5),
            'alpha':round(al,3),'beta':round(be,3),'gamma':round(ga,3),
            'system':            (ph.get('system') or 'triclinic').lower(),
            'spacegroup_number': ph.get('spacegroup_number', 1),
            'spacegroup':        ph.get('spacegroup', ''),
            'scale':             round(float(st['S']), 5),
            'U':round(st['U'],5),'V':round(st['V'],5),
            'W':round(st['W'],5),
            'X':round(st['X'],5),'Y':round(st['Y'],5),
            'eta_at_strongest':  round(eta_sc, 3),
            'fwhm_deg':          round(fwhm_sc, 4),
            'crystallite_size_A':  round(cryst_A, 1) if cryst_A else None,
            'crystallite_size_nm': round(cryst_A/10, 2) if cryst_A else None,
            'weight_fraction_%':   round(wt_frac, 1),
            'n_reflections':       len(st['refs']),
            'tick_positions':      _filter_tick_positions(st['refs'], st['I_hkl']),
            'seeded_by':           st.get('seeded', 'unknown'),
        })

    return {
        'tt':             tt_r.tolist(),
        'y_obs':          y_r.tolist(),
        'y_calc':         yc_f.tolist(),
        'y_background':   bg_f.tolist(),
        'phase_patterns': phase_patterns,
        'residuals':      diff_f.tolist(),
        'statistics':     stats,
        'phase_results':  phase_results,
        'zero_shift':     round(float(zero), 5),
        'n_params':       n_params,
        'wavelength':     wavelength,
        'global_scale':   round(float(sum(st['S'] for st in phase_state) / max(len(phase_state),1)), 5),
        'pymatgen_used':  any(st.get('seeded') == 'pymatgen' for st in phase_state),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# RIETVELD REFINEMENT
# ═══════════════════════════════════════════════════════════════════════════════

def run_rietveld(tt, y_obs, sigma, phases, wavelength,
                 tt_min=None, tt_max=None,
                 n_bg_coeffs=6, max_iter=50,
                 progress_callback=None):
    """
    Rietveld refinement — structure-constrained profile fitting.

    Unlike Le Bail (where each peak intensity is free), Rietveld computes
    peak intensities from the crystal structure:
        I_hkl = S × mult × |F(hkl)|² × LP × DW(B_iso)

    This dramatically reduces the parameter count and constrains the
    relative intensities, fixing the high-angle underfitting problem.

    Refined parameters per phase: S (scale), a [b,c...] (cell), U,V,W,η (profile), B_iso
    Global parameters: background (Chebyshev), zero shift

    Parameters
    ----------
    tt, y_obs, sigma : np.ndarray — data arrays
    phases   : list of dicts — must include 'cif_text' or 'sites' for atom positions
    wavelength : float (Å)
    tt_min/max : float
    n_bg_coeffs : int
    max_iter    : int

    Returns dict compatible with Le Bail output (same keys for plotting).
    """
    if tt_min is None: tt_min = float(tt.min())
    if tt_max is None: tt_max = float(tt.max())

    mask   = (tt >= tt_min) & (tt <= tt_max)
    tt_r   = tt[mask]
    y_r    = y_obs[mask]
    sig_r  = sigma[mask] if sigma is not None else np.sqrt(np.maximum(y_r, 1.0))
    weights = 1.0 / np.maximum(sig_r**2, 1e-6)
    n_pts  = len(tt_r)

    if progress_callback:
        progress_callback('Rietveld: generating reflections...')

    # ── Build per-phase state ─────────────────────────────────────────────
    bg_est = float(np.percentile(y_r, 5))
    phase_state = []

    for ph in phases:
        # Extract atom sites — use pymatgen expansion for correct F²,
        # then fall back to raw parse_cif (asymmetric unit only).
        sites = ph.get('sites')
        if not sites and ph.get('cif_text'):
            sites = expand_sites_from_cif(ph['cif_text'])
            if not sites:
                try:
                    parsed = _parse_cif_cryst(ph['cif_text'])
                    sites = parsed.get('sites') or None
                except Exception:
                    sites = None
        if not sites:
            raise ValueError(f"Phase '{ph.get('name','?')}' has no atom sites. "
                             f"Rietveld requires a CIF with atomic coordinates.")

        free_v, free_n = _cell_free(ph)
        a, b, c, al, be, ga = _full_cell(free_v, free_n, ph)
        sys_ = (ph.get('system') or 'triclinic').lower()
        sg   = ph.get('spacegroup_number', 1)

        refs = generate_reflections_rietveld(
            a, b, c, al, be, ga, sys_, sg,
            wavelength, tt_min, tt_max, sites, hkl_max=12
        )

        # Initial B_iso estimate (0.5 Å² is typical for metals)
        B_init = ph.get('B_iso_init', 0.5)

        phase_state.append({
            'ph':     ph,
            'sites':  sites,
            'refs':   refs,
            'S':      1.0,
            'B_iso':  B_init,
            'U':      ph.get('U_init', 0.01),
            'V':      ph.get('V_init', -0.01),
            'W':      ph.get('W_init', 0.15),
            'X':      ph.get('X_init', 0.0),
            'Y':      ph.get('Y_init', 0.1),
            'eta':    ph.get('eta_init', 0.5),
            'free_v': free_v,
            'free_n': free_n,
            'sys':    sys_,
            'sg':     sg,
        })

    # ── Background and zero shift ─────────────────────────────────────────
    n_bg  = n_bg_coeffs
    bg_c  = np.zeros(n_bg)
    bg_c[0] = bg_est
    zero  = 0.0

    # ── Helper: compute full y_calc from current parameters ───────────────
    def calc_pattern(states, bg_c_, zero_):
        bg = np.maximum(chebyshev_background(tt_r, bg_c_, tt_min, tt_max), 0)
        pat = np.zeros(n_pts)
        for st in states:
            I_hkl = compute_rietveld_intensities(
                st['refs'], st['sites'], {'_all': st['B_iso']})
            profs = _get_profiles(tt_r, _refs_to_legacy(st['refs']),
                                   st['U'], st['V'], st['W'],
                                   st.get('eta', 0.5), zero_,
                                   st['X'], st['Y'])
            for k in range(len(st['refs'])):
                pat += st['S'] * I_hkl[k] * profs[k]
        return pat + bg, bg, pat

    def rwp(yo, yc, w):
        d = yo - yc
        return math.sqrt(np.sum(w*d**2) / np.sum(w*yo**2)) * 100

    # ── Better initial scale estimate ───────────────────────────────────
    # Compute theoretical pattern at B_iso=0.5 and match to observed data
    for st in phase_state:
        I_hkl = compute_rietveld_intensities(
            st['refs'], st['sites'], {'_all': st['B_iso']})
        profs = _get_profiles(tt_r, _refs_to_legacy(st['refs']),
                               st['U'], st['V'], st['W'],
                               st.get('eta', 0.5), zero,
                               st['X'], st['Y'])
        this_pat = np.zeros(n_pts)
        for k in range(len(st['refs'])):
            this_pat += I_hkl[k] * profs[k]
        bg_init = np.maximum(chebyshev_background(tt_r, bg_c, tt_min, tt_max), 0)
        resid_target = np.maximum(y_r - bg_init, 0.1)
        S_num = np.sum(weights * resid_target * this_pat)
        S_den = np.sum(weights * this_pat**2)
        st['S'] = max(S_num / max(S_den, 1e-10), 1e-6)

    if progress_callback:
        progress_callback('Rietveld: refining...')

    # ── Main iteration loop with staged parameter release ─────────────
    # Stage 1 (iter 0-2):  scale + profile + zero  (B_iso & cell FIXED)
    # Stage 2 (iter 3-6):  + cell parameters
    # Stage 3 (iter 7+):   + B_iso  (all free)
    prev_rwp = 999.0

    for iteration in range(max_iter):

        # Determine which parameters are active this iteration
        refine_cell = (iteration >= 3)
        refine_biso = (iteration >= 7)

        # ── Step 1: Update per-phase scale factors (linear) ───────────────
        bg_cur = np.maximum(chebyshev_background(tt_r, bg_c, tt_min, tt_max), 0)
        for i_ph, st in enumerate(phase_state):
            I_hkl = compute_rietveld_intensities(
                st['refs'], st['sites'], {'_all': st['B_iso']})
            profs = _get_profiles(tt_r, _refs_to_legacy(st['refs']),
                                   st['U'], st['V'], st['W'],
                                   st.get('eta', 0.5), zero,
                                   st['X'], st['Y'])
            other = np.zeros(n_pts)
            for j, st2 in enumerate(phase_state):
                if j == i_ph:
                    continue
                I2 = compute_rietveld_intensities(
                    st2['refs'], st2['sites'], {'_all': st2['B_iso']})
                p2 = _get_profiles(tt_r, _refs_to_legacy(st2['refs']),
                                    st2['U'], st2['V'], st2['W'],
                                    st2.get('eta', 0.5), zero,
                                    st2['X'], st2['Y'])
                for k in range(len(st2['refs'])):
                    other += st2['S'] * I2[k] * p2[k]

            this_pat = np.zeros(n_pts)
            for k in range(len(st['refs'])):
                this_pat += I_hkl[k] * profs[k]
            resid_target = y_r - bg_cur - other
            S_num = np.sum(weights * resid_target * this_pat)
            S_den = np.sum(weights * this_pat**2)
            st['S'] = max(S_num / max(S_den, 1e-10), 1e-6)

        # ── Step 2: Update background ─────────────────────────────────────
        yc_nobg, _, pat_total = calc_pattern(phase_state, np.zeros(n_bg), zero)
        pat_total_only = pat_total
        def resid_bg(x):
            bg_ = np.maximum(chebyshev_background(tt_r, x, tt_min, tt_max), 0)
            return (y_r - pat_total_only - bg_) * np.sqrt(weights)
        rb = least_squares(resid_bg, bg_c, method='trf',
                           max_nfev=200, ftol=1e-10, verbose=0)
        bg_c = rb.x

        # ── Step 3: SIMULTANEOUS multi-phase refinement (staged) ──────────
        # All phases' parameters are refined together so that overlapping
        # peaks are correctly apportioned between phases.
        bg_fixed = np.maximum(chebyshev_background(tt_r, bg_c, tt_min, tt_max), 0)

        # Build layout: per-phase parameter offsets
        # Each phase block: [cell_params...?, S, B_iso?, U, V, W, X, Y]
        # Followed by global [zero]
        phase_layouts = []  # list of (offset, n_cell, has_biso) per phase
        x0, lo, hi = [], [], []
        _ref_caches = [{'fv': None, 'refs': st['refs']} for st in phase_state]

        for i_ph, st in enumerate(phase_state):
            ph = st['ph']
            free_n = st['free_n']
            offset = len(x0)

            n_cell = 0
            if refine_cell:
                fv0 = list(st['free_v'])
                n_cell = len(fv0)
                x0.extend(fv0)
                lo.extend([v * 0.92 for v in fv0])
                hi.extend([v * 1.08 for v in fv0])
                for i_p, nm in enumerate(free_n):
                    if nm in ('alpha', 'beta', 'gamma'):
                        lo[offset + i_p] = fv0[i_p] - 5
                        hi[offset + i_p] = fv0[i_p] + 5

            # Scale
            x0.append(st['S']); lo.append(1e-6); hi.append(1e4)

            # B_iso
            has_biso = refine_biso
            if refine_biso:
                x0.append(st['B_iso']); lo.append(0.0); hi.append(8.0)

            # U, V, W, X, Y
            x0.extend([st['U'], st['V'], st['W'], st['X'], st['Y']])
            lo.extend([0.0, -1.0, 0.005, 0.0, 0.0])
            hi.extend([5.0,  1.0, 3.0,   2.0, 5.0])

            phase_layouts.append((offset, n_cell, has_biso))

        # Global zero
        x0.append(zero); lo.append(-0.5); hi.append(0.5)

        def resid_joint(x):
            zero_ = x[-1]
            pat = np.zeros(n_pts)

            for i_ph, st in enumerate(phase_state):
                off, n_cell, has_biso = phase_layouts[i_ph]
                idx = off
                ph = st['ph']
                free_n = st['free_n']

                if refine_cell and n_cell > 0:
                    fv = x[idx:idx + n_cell]
                    idx += n_cell
                else:
                    fv = list(st['free_v'])

                S_ph = max(x[idx], 1e-6); idx += 1
                B_iso = max(x[idx], 0.0) if has_biso else st['B_iso']
                if has_biso: idx += 1

                U_ = x[idx]; V_ = x[idx+1]; W_ = max(x[idx+2], 0.005)
                X_ = x[idx+3]; Y_ = x[idx+4]

                cache = _ref_caches[i_ph]
                if refine_cell and n_cell > 0:
                    if (cache['fv'] is None or
                            np.max(np.abs(np.array(fv) - np.array(cache['fv']))) > 1e-4):
                        a_,b_,c_,al_,be_,ga_ = _full_cell(fv, free_n, ph)
                        try:
                            cache['refs'] = generate_reflections_rietveld(
                                a_,b_,c_,al_,be_,ga_, st['sys'], st['sg'],
                                wavelength, tt_min, tt_max, st['sites'], hkl_max=12)
                        except Exception:
                            pass
                        cache['fv'] = list(fv)

                refs_l = cache['refs']
                I_hkl = compute_rietveld_intensities(refs_l, st['sites'], {'_all': B_iso})
                profs = _get_profiles(tt_r, _refs_to_legacy(refs_l),
                                       U_, V_, W_, 0.5, zero_, X_, Y_)
                for k in range(len(refs_l)):
                    pat += S_ph * I_hkl[k] * profs[k]

            return (y_r - pat - bg_fixed) * np.sqrt(weights)

        try:
            rs = least_squares(resid_joint, x0, bounds=(lo, hi),
                               method='trf', max_nfev=800,
                               ftol=1e-7, xtol=1e-7, verbose=0)
            x_out = rs.x
        except Exception:
            x_out = np.array(x0)

        # Unpack results back into phase_state
        zero = x_out[-1]
        for i_ph, st in enumerate(phase_state):
            off, n_cell, has_biso = phase_layouts[i_ph]
            ph = st['ph']
            free_n = st['free_n']
            sys_ = st['sys']
            sg = st['sg']
            sites = st['sites']
            idx = off

            if refine_cell and n_cell > 0:
                fv_new = x_out[idx:idx + n_cell]
                idx += n_cell
                st['free_v'] = fv_new
                for nm, vl in zip(free_n, fv_new):
                    ph[nm] = vl
                a_n, b_n, c_n, al_n, be_n, ga_n = _full_cell(fv_new, free_n, ph)
                if sys_ == 'cubic':
                    ph['b'] = ph['c'] = a_n
                elif sys_ in ('hexagonal', 'trigonal', 'tetragonal'):
                    ph['b'] = a_n
                ph['c'] = c_n
                try:
                    st['refs'] = generate_reflections_rietveld(
                        a_n, b_n, c_n, al_n, be_n, ga_n, sys_, sg,
                        wavelength, tt_min, tt_max, sites, hkl_max=12)
                except Exception:
                    pass

            st['S'] = max(x_out[idx], 1e-6); idx += 1
            if has_biso:
                st['B_iso'] = max(x_out[idx], 0.0); idx += 1
            st['U'] = x_out[idx]; st['V'] = x_out[idx+1]
            st['W'] = max(x_out[idx+2], 0.005)
            st['X'] = x_out[idx+3]; st['Y'] = x_out[idx+4]

        # ── Convergence check ─────────────────────────────────────────────
        yc_cur, _, _ = calc_pattern(phase_state, bg_c, zero)
        cur_rwp = rwp(y_r, yc_cur, weights)
        stage_name = 'profile' if not refine_cell else ('cell' if not refine_biso else 'all')
        if progress_callback:
            progress_callback(
                f'Rietveld iter {iteration+1} [{stage_name}]: Rwp = {cur_rwp:.2f}%')
        # Only check convergence after all parameters are free (stage 3)
        if refine_biso and abs(prev_rwp - cur_rwp) < 0.1:
            break
        prev_rwp = cur_rwp

    # ── Final results (same format as Le Bail for compatibility) ──────────
    if progress_callback:
        progress_callback('Rietveld: computing final statistics...')

    yc_f, bg_f, pat_f = calc_pattern(phase_state, bg_c, zero)
    diff_f = y_r - yc_f
    n_params = (sum(len(st['free_n']) + 7 for st in phase_state)  # cell+S+B+U+V+W+X+Y
                + n_bg + 1)  # bg + zero
    stats = compute_fit_statistics(y_r, yc_f, weights, n_params)

    phase_patterns = []
    phase_results  = []

    # ── Weight fractions via Hill & Howard (1987) ────────────────────────
    zmv_values_r = []
    use_zmv_r = True
    for st in phase_state:
        ph = st['ph']
        a_ = float(ph.get('a', 4.0))
        b_ = float(ph.get('b') or a_)
        c_ = float(ph.get('c') or a_)
        al_ = float(ph.get('alpha', 90.0) or 90.0)
        be_ = float(ph.get('beta',  90.0) or 90.0)
        ga_ = float(ph.get('gamma', 90.0) or 90.0)
        V = cell_volume(a_, b_, c_, al_, be_, ga_)
        Z = ph.get('Z')
        M = molar_mass_from_formula(ph.get('formula', ''))
        if Z and M:
            zmv_values_r.append(float(st['S']) * float(Z) * float(M) * V)
        else:
            use_zmv_r = False
            break

    if not use_zmv_r:
        zmv_values_r = []
        for st in phase_state:
            I_hkl_tmp = compute_rietveld_intensities(
                st['refs'], st['sites'], {'_all': st['B_iso']})
            profs_tmp = _get_profiles(tt_r, _refs_to_legacy(st['refs']),
                                       st['U'], st['V'], st['W'],
                                       st.get('eta', 0.5), zero,
                                       st['X'], st['Y'])
            area = st['S'] * sum(
                I_hkl_tmp[k] * profs_tmp[k].sum() for k in range(len(st['refs'])))
            zmv_values_r.append(float(area))

    total_zmv_r = sum(zmv_values_r) or 1e-10

    for i_ph, st in enumerate(phase_state):
        I_hkl = compute_rietveld_intensities(
            st['refs'], st['sites'], {'_all': st['B_iso']})
        # Use tight profiles (3× FWHM) for display to prevent Lorentzian
        # tails from one phase bleeding into another phase's peak regions.
        display_profs = _get_profiles(tt_r, _refs_to_legacy(st['refs']),
                                       st['U'], st['V'], st['W'],
                                       st.get('eta', 0.5), zero,
                                       st['X'], st['Y'],
                                       window_factor=3.0)
        pat_ph = st['S'] * sum(I_hkl[k] * display_profs[k]
                               for k in range(len(st['refs'])))
        phase_patterns.append(pat_ph.tolist())

        ph = st['ph']
        a  = float(ph.get('a', 4.0))
        b  = float(ph.get('b') or a)
        c  = float(ph.get('c') or a)
        al = float(ph.get('alpha', 90.0) or 90.0)
        be = float(ph.get('beta',  90.0) or 90.0)
        ga = float(ph.get('gamma', 90.0) or 90.0)

        strongest = max(st['refs'], key=lambda r: r['mult']) if st['refs'] else None
        strongest_tt = strongest['two_theta'] if strongest else 39.0

        # Crystallite size: prefer TCH Y parameter (direct size extraction)
        Y_val = st['Y']
        if Y_val > 0.01:
            cryst_A = size_from_Y(Y_val, wavelength)
        else:
            fwhm_tot, _ = tch_fwhm_eta(strongest_tt, st['U'], st['V'],
                                         st['W'], st['X'], st['Y'])
            cryst_A = scherrer_size(fwhm_tot, strongest_tt, wavelength)

        fwhm_sc, eta_sc = tch_fwhm_eta(strongest_tt, st['U'], st['V'],
                                         st['W'], st['X'], st['Y'])

        wt_frac = (zmv_values_r[i_ph] / total_zmv_r) * 100

        _sg_sym_rv  = (ph.get('spacegroup', '') or
                       f"SG{ph.get('spacegroup_number', 1)}")
        _base_rv    = ph.get('name') or ph.get('formula') or f'Phase {i_ph+1}'
        if _sg_sym_rv and _sg_sym_rv.replace(' ','') not in _base_rv.replace(' ',''):
            _display_rv = f"{_base_rv} {_sg_sym_rv}"
        else:
            _display_rv = _base_rv

        phase_results.append({
            'name':              _display_rv,
            'cod_id':            ph.get('cod_id', ''),
            'formula':           ph.get('formula', ''),
            'a': round(a,5), 'b': round(b,5), 'c': round(c,5),
            'alpha': round(al,3), 'beta': round(be,3), 'gamma': round(ga,3),
            'system':            (ph.get('system') or 'triclinic').lower(),
            'spacegroup_number': ph.get('spacegroup_number', 1),
            'spacegroup':        ph.get('spacegroup', ''),
            'scale':             round(float(st['S']), 5),
            'B_iso':             round(float(st['B_iso']), 4),
            'U': round(st['U'],5), 'V': round(st['V'],5),
            'W': round(st['W'],5),
            'X': round(st['X'],5), 'Y': round(st['Y'],5),
            'eta_at_strongest':  round(eta_sc, 3),
            'fwhm_deg':          round(fwhm_sc, 4),
            'crystallite_size_A':  round(cryst_A, 1) if cryst_A else None,
            'crystallite_size_nm': round(cryst_A/10, 2) if cryst_A else None,
            'weight_fraction_%':   round(wt_frac, 1),
            'n_reflections':       len(st['refs']),
            'tick_positions':      _filter_tick_positions(st['refs'], I_hkl),
            'seeded_by':           'rietveld',
        })

    return {
        'tt':             tt_r.tolist(),
        'y_obs':          y_r.tolist(),
        'y_calc':         yc_f.tolist(),
        'y_background':   bg_f.tolist(),
        'phase_patterns': phase_patterns,
        'residuals':      diff_f.tolist(),
        'statistics':     stats,
        'phase_results':  phase_results,
        'zero_shift':     round(float(zero), 5),
        'n_params':       n_params,
        'wavelength':     wavelength,
        'global_scale':   round(float(sum(st['S'] for st in phase_state) / max(len(phase_state),1)), 5),
        'pymatgen_used':  False,
        'method':         'rietveld',
    }


def _refs_to_legacy(riet_refs):
    """
    Convert Rietveld-style ref dicts to legacy tuples for _get_profiles().
    _get_profiles expects [(two_theta, d, hkl, mult), ...]
    """
    return [(r['two_theta'], r['d'], r['hkl'], r['mult']) for r in riet_refs]
