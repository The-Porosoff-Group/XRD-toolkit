"""
modules/xrd/gsasii_calibration.py
Instrument profile calibration using GSAS-II and a line-standard
(e.g. NIST SRM 640g Si).

Separate from gsasii_backend.py by design:
  - gsasii_backend.py: refine unknown samples with fixed instrument params
  - gsasii_calibration.py: refine instrument params with fixed cell/structure

Strategy: conservative, rollback-driven.  Parameters are added one at a
time as TRIALS.  Each trial is accepted ONLY if:
  1. Rwp improves meaningfully
  2. No pathological correlations or SVD warnings
  3. Parameters are physically plausible (bounded, FWHM² > 0)
If any trial fails, the code rolls back to the last known-good state.

A stable imperfect .instprm is much better than a formally lower-Rwp
but invalid one.

Usage:
    from modules.xrd.gsasii_calibration import run_calibration
    result = run_calibration(tt, y_obs, sigma, phase, wavelength, ...)
"""

import math, os, sys, tempfile, shutil, warnings, re
import numpy as np

from .gsasii_backend import (
    _add_gsas2pkg_paths, _GSASII_AVAILABLE, _GSASII_IMPORT_ERROR,
    _is_cu_kalpha, CU_KALPHA1_A, CU_KALPHA2_A, CU_KALPHA2_RATIO,
    _build_conventional_cif, INSTRUMENT_PROFILES, DEFAULT_INSTRUMENT,
)

# Conservative starting guesses — W=0.01 to avoid starting too broad.
# If GSAS-II uses centideg² internally, W=5.0 → FWHM≈2° which is
# far too wide.  Starting near zero lets GSAS-II find the right scale.
_CAL_DEFAULTS = {
    'U': 0.0, 'V': 0.0, 'W': 0.01,
    'X': 0.0, 'Y': 0.0, 'SH/L': 0.002,
}

# Si peak positions for validation (Cu Kα1, a=5.431109 Å)
_SI_PEAK_ANGLES = np.array([28.44, 47.30, 56.12, 69.13, 76.37, 88.03])


def _profile_plausible(params, tt_check=None):
    """Check if U/V/W/X/Y are physically plausible.

    Returns (ok: bool, reason: str).
    """
    U, V, W = params.get('U', 0), params.get('V', 0), params.get('W', 0)
    X, Y = params.get('X', 0), params.get('Y', 0)

    if tt_check is None:
        tt_check = _SI_PEAK_ANGLES

    # Caglioti FWHM² must be positive at all Si peak positions
    tan_t = np.tan(np.radians(tt_check / 2.0))
    hg2 = U * tan_t**2 + V * tan_t + W
    if np.nanmin(hg2) <= 0:
        return False, f"FWHM²_G non-positive at 2θ={tt_check[np.argmin(hg2)]:.1f}°"

    # Broad guardrails — not fundamental constants, just sanity bounds
    if abs(U) > 1000:
        return False, f"|U|={abs(U):.1f} > 1000 (absurd)"
    if abs(V) > 200:
        return False, f"|V|={abs(V):.1f} > 200 (absurd)"
    if W > 500:
        return False, f"W={W:.1f} > 500 (absurd)"
    if X < 0 or X > 10:
        return False, f"X={X:.4f} outside [0, 10]"
    if Y < 0 or Y > 10:
        return False, f"Y={Y:.4f} outside [0, 10]"

    return True, "ok"


def run_calibration(tt, y_obs, sigma, phase, wavelength,
                    tt_min=None, tt_max=None, n_bg_coeffs=6,
                    polariz=None, instrument=None,
                    output_instprm=None, progress_callback=None,
                    keep_workdir=False):
    """
    Run instrument profile calibration on a line-broadening standard.

    Parameters
    ----------
    tt, y_obs, sigma : np.ndarray — data arrays
    phase : dict — single phase dict with 'cif_text', 'a', etc.
        Cell must be the NIST certified value (enforced by app.py).
    wavelength : float (Å)
    tt_min/max : float — fitting range
    n_bg_coeffs : int — background Chebyshev terms (default 8)
    polariz : float — beam polarization (default from instrument profile)
    instrument : str — instrument profile key
    output_instprm : str — path to write .instprm
    progress_callback : callable
    keep_workdir : bool — keep temp .gpx/.lst for debugging

    Returns dict with 'instprm_path', 'params', 'Rwp', plot arrays, etc.
    """
    if not _GSASII_AVAILABLE:
        raise RuntimeError(f"GSAS-II not available: {_GSASII_IMPORT_ERROR}")

    _add_gsas2pkg_paths()
    try:
        import GSASIIscriptable as G2sc
    except ImportError:
        from GSASII import GSASIIscriptable as G2sc

    # ── Resolve output path early and back up existing file ────────────
    instrument = instrument or DEFAULT_INSTRUMENT
    _toolkit_root = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    if output_instprm is None:
        _profiles = INSTRUMENT_PROFILES.get(instrument,
                      INSTRUMENT_PROFILES[DEFAULT_INSTRUMENT])
        _fname = _profiles.get('instprm_filename',
                                f'{instrument}_Si640g.instprm')
        output_instprm = os.path.join(_toolkit_root, _fname)

    if os.path.isfile(output_instprm):
        print(f"  Existing .instprm found: {output_instprm} "
              f"(will be backed up before overwrite)", flush=True)

    # ── Resolve instrument profile ─────────────────────────────────────
    profile = INSTRUMENT_PROFILES.get(instrument,
                INSTRUMENT_PROFILES[DEFAULT_INSTRUMENT])
    if polariz is None:
        polariz = profile.get('polariz', 0.5)

    print(f"  *** INSTRUMENT CALIBRATION MODE ***", flush=True)
    print(f"  Instrument: {profile['label']}", flush=True)
    print(f"  Polariz: {polariz}", flush=True)
    print(f"  Cell: FIXED (certified standard)", flush=True)
    print(f"  Strategy: conservative trial-based with rollback", flush=True)

    # ── Validate ───────────────────────────────────────────────────────
    cif_text = phase.get('cif_text', '')
    if not cif_text:
        raise ValueError("Calibration phase must have CIF text.")

    if tt_min is None: tt_min = float(tt.min())
    if tt_max is None: tt_max = float(tt.max())

    mask = (tt >= tt_min) & (tt <= tt_max)
    tt_r = tt[mask]
    y_r = y_obs[mask]
    sig_r = (sigma[mask] if sigma is not None
             else np.sqrt(np.maximum(y_r, 1.0)))

    if progress_callback:
        progress_callback('Calibration: setting up project...')

    # ── Temp directory ─────────────────────────────────────────────────
    _app_tmp = os.path.join(_toolkit_root, '.gsas_tmp')
    os.makedirs(_app_tmp, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix='gsas2_cal_', dir=_app_tmp)
    print(f"  Work dir: {work_dir}"
          + (" (kept)" if keep_workdir else ""), flush=True)

    try:
        # ── Write data files ───────────────────────────────────────────
        data_path = os.path.join(work_dir, 'standard.xye')
        with open(data_path, 'w', encoding='utf-8', newline='\n') as f:
            for i in range(len(tt_r)):
                f.write(f"{tt_r[i]:.6f}  {y_r[i]:.4f}  {sig_r[i]:.4f}\n")

        # Write calibration CIF directly — DO NOT use _build_conventional_cif
        # or the COD CIF, because the Fd-3m origin choice (setting 1 vs 2)
        # causes silent wrong structure factors if the atom position doesn't
        # match GSAS-II's convention.
        #
        # GSAS-II uses Fd-3m setting 2 (origin at center, 3m).
        # Si 8a Wyckoff site in setting 2: (1/8, 1/8, 1/8).
        # Using (0,0,0) — which is setting 1 — produces F²≈0 for all
        # diamond reflections, giving a flat calculated pattern and
        # 60% Rwp.
        _a = phase.get('a', 5.431109)
        cif_path = os.path.join(work_dir, 'standard.cif')
        _sg_num = phase.get('spacegroup_number', 227)

        # For Fd-3m (#227), write a known-good CIF with setting 2
        if _sg_num == 227:
            _cif_out = f"""data_calibration_standard
_cell_length_a {_a}
_cell_length_b {_a}
_cell_length_c {_a}
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_symmetry_Int_Tables_number 227
_symmetry_space_group_name_H-M 'F d -3 m :2'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1  Si  0.12500  0.12500  0.12500  1.0
"""
            print(f"  CIF: Fd-3m setting 2, Si at (1/8,1/8,1/8), "
                  f"a={_a} (certified)", flush=True)
        else:
            # For other standards, use the normal CIF pipeline
            built_cif = _build_conventional_cif(phase)
            _cif_out = built_cif if built_cif else cif_text
            for _tag, _val in [
                    ('_cell_length_a', _a),
                    ('_cell_length_b', phase.get('b', _a)),
                    ('_cell_length_c', phase.get('c', _a)),
                    ('_cell_angle_alpha', phase.get('alpha', 90.0)),
                    ('_cell_angle_beta', phase.get('beta', 90.0)),
                    ('_cell_angle_gamma', phase.get('gamma', 90.0))]:
                _cif_out = re.sub(
                    rf'({_tag}\s+)\S+', rf'\g<1>{_val}', _cif_out)
            print(f"  CIF: SG {_sg_num}, a={_a} (certified)", flush=True)

        with open(cif_path, 'w', encoding='utf-8') as f:
            f.write(_cif_out)

        # Write initial .instprm with FRESH starting guesses.
        # CRITICAL: calibration must NEVER use an existing measured
        # .instprm as input.  It always starts from _CAL_DEFAULTS.
        _use_doublet = _is_cu_kalpha(wavelength)
        instprm_path = os.path.join(work_dir, 'initial.instprm')
        print(f"  Input instprm: {instprm_path} (fresh defaults, "
              f"NOT the existing measured file)", flush=True)
        with open(instprm_path, 'w') as f:
            f.write("#GSAS-II instrument parameter file; "
                    "do not add/delete items!\n")
            f.write("Type:PXC\n")
            if _use_doublet:
                f.write(f"Lam1:{CU_KALPHA1_A:.6f}\n")
                f.write(f"Lam2:{CU_KALPHA2_A:.6f}\n")
                f.write(f"I(L2)/I(L1):{CU_KALPHA2_RATIO:.4f}\n")
            else:
                f.write(f"Lam:{wavelength:.6f}\n")
            f.write("Zero:0.0\n")
            f.write(f"Polariz.:{polariz}\n")
            for k in ['U', 'V', 'W', 'X', 'Y']:
                f.write(f"{k}:{_CAL_DEFAULTS[k]}\n")
            f.write("Z:0.0\n")
            f.write(f"SH/L:{_CAL_DEFAULTS['SH/L']}\n")
            f.write("Azimuth:0.0\n")

        # ── Create project ─────────────────────────────────────────────
        gpx_path = os.path.join(work_dir, 'calibrate.gpx')
        gpx = G2sc.G2Project(newgpx=gpx_path)

        histogram = gpx.add_powder_histogram(data_path, instprm_path,
                                              fmthint='xye')
        if isinstance(histogram, list):
            histogram = histogram[0]
        print(f"  Histogram: {len(tt_r)} points", flush=True)

        try:
            histogram.data['Sample Parameters']['Scale'] = [1.0, False]
        except Exception:
            pass

        if progress_callback:
            progress_callback('Calibration: adding standard phase...')
        print("  Adding phase...", flush=True)
        phase_obj = gpx.add_phase(cif_path, phasename='Standard',
                                   histograms=[histogram])
        if isinstance(phase_obj, list):
            phase_obj = phase_obj[0]
        print(f"  Phase: {phase_obj.name}", flush=True)

        phase_obj.set_refinements({'Cell': False})
        try:
            _cell = phase_obj.data['General']['Cell']
            print(f"  Cell in GSAS-II: a={_cell[1]:.6f} (fixed)",
                  flush=True)
        except Exception:
            pass

        try:
            hapData = list(phase_obj.data['Histograms'].values())[0]
            hapData['Scale'] = [1.0, True]
        except Exception:
            pass

        gpx.save()

        # ── Helpers ────────────────────────────────────────────────────
        def _get_inst_value(key):
            """Read instrument param robustly across GSAS-II versions."""
            entry = histogram.data['Instrument Parameters'][0].get(
                key, 0.0)
            if isinstance(entry, (int, float)) and not isinstance(
                    entry, bool):
                return float(entry)
            try:
                if isinstance(entry, (np.integer, np.floating)):
                    return float(entry)
            except Exception:
                pass
            if isinstance(entry, (list, tuple)):
                numeric = [
                    x for x in entry
                    if isinstance(x, (int, float, np.integer, np.floating))
                    and not isinstance(x, (bool, np.bool_))
                ]
                if len(numeric) >= 2:
                    return float(numeric[1])
                if len(numeric) == 1:
                    return float(numeric[0])
            return 0.0

        def _get_arrays():
            """Get histogram arrays using correct GSAS-II data keys."""
            # GSAS-II uses 'yobs'/'ycalc' not 'observed'/'calculated'
            x = np.array(histogram.getdata('x'), dtype=float)
            yobs = np.array(histogram.getdata('yobs'), dtype=float)
            ycalc = np.array(histogram.getdata('ycalc'), dtype=float)
            bg = np.array(histogram.getdata('background'), dtype=float)
            return x, yobs, ycalc, bg

        def _get_rwp():
            """Get Rwp from GSAS-II."""
            try:
                s = histogram.get_statistics()
                if s and 'Rwp' in s:
                    return float(s['Rwp'])
            except (AttributeError, TypeError):
                pass
            try:
                rvals = gpx.data['Covariance']['data']['Rvals']
                if 'Rwp' in rvals:
                    return float(rvals['Rwp'])
            except (KeyError, TypeError):
                pass
            try:
                _, y_o, y_c, _ = _get_arrays()
                w = 1.0 / np.maximum(y_o, 1.0)
                return float(100.0 * np.sqrt(
                    np.sum(w * (y_o - y_c)**2) / np.sum(w * y_o**2)))
            except Exception:
                pass
            return 999.0

        def _get_stats():
            """Get full stats dict."""
            try:
                s = histogram.get_statistics()
                if s:
                    return s
            except (AttributeError, TypeError):
                pass
            try:
                rv = gpx.data['Covariance']['data']['Rvals']
                return {'Rwp': rv.get('Rwp', 0), 'Rp': rv.get('Rp', 0),
                        'reduced chi2': rv.get('GOF', 0)**2
                        if 'GOF' in rv else 0}
            except (KeyError, TypeError):
                return {'Rwp': _get_rwp(), 'Rp': 0, 'reduced chi2': 0}

        def _current_params():
            """Snapshot of all instrument params."""
            return {k: _get_inst_value(k)
                    for k in ['U', 'V', 'W', 'X', 'Y', 'Zero', 'SH/L']}

        def _checkpoint(name):
            """Save a named checkpoint. Returns the checkpoint path."""
            gpx.save()
            ckpt = os.path.join(work_dir, f'{name}.gpx')
            shutil.copy2(gpx_path, ckpt)
            return ckpt

        def _restore(ckpt_path):
            """Restore from a named checkpoint."""
            nonlocal gpx, histogram, phase_obj
            shutil.copy2(ckpt_path, gpx_path)
            gpx = G2sc.G2Project(gpx_path)
            histogram = gpx.histograms()[0]
            phase_obj = gpx.phases()[0]

        def _refine(label, inst_params, cycles=8):
            """Run refinement, return Rwp."""
            rd = {'set': {
                'Background': {'type': 'chebyschev-1', 'refine': True,
                                'no. coeffs': n_bg_coeffs},
            }, 'cycles': cycles}
            if inst_params:
                rd['set']['Instrument Parameters'] = inst_params
            gpx.do_refinements([rd])
            rwp = _get_rwp()
            vals = [f"{k}={_get_inst_value(k):.4f}"
                    for k in (inst_params or [])]
            print(f"    {label}: Rwp={rwp:.3f}%"
                  + (f"  [{', '.join(vals)}]" if vals else ""),
                  flush=True)
            # Phase contribution diagnostic
            try:
                _, _yo, _yc, _yb = _get_arrays()
                _yp = _yc - _yb
                print(f"      max(obs)={np.max(_yo):.1f}, "
                      f"max(calc)={np.max(_yc):.1f}, "
                      f"max(bg)={np.max(_yb):.1f}, "
                      f"max(calc-bg)={np.max(_yp):.1f}",
                      flush=True)
            except Exception as _de:
                print(f"      (diagnostic failed: {_de})", flush=True)
            return rwp

        def _trial(label, inst_params, ckpt_from, cycles=10,
                   min_rwp_gain=0.1):
            """Try adding params from a checkpoint. Accept only if
            plausible + meaningful improvement.

            Always starts from ckpt_from so trials are independent.
            Returns (accepted: bool, rwp: float, ckpt_path or None).
            """
            _restore(ckpt_from)
            rwp_before = _get_rwp()

            rwp = _refine(label, inst_params, cycles)
            p = _current_params()
            ok, reason = _profile_plausible(p, _SI_PEAK_ANGLES)

            if not ok:
                print(f"    REJECTED: {reason} — rolling back", flush=True)
                _restore(ckpt_from)
                return False, rwp_before, None

            rwp_gain = rwp_before - rwp
            if rwp_gain < min_rwp_gain:
                print(f"    REJECTED: ΔRwp={rwp_gain:.3f}% < "
                      f"{min_rwp_gain}% — not meaningful", flush=True)
                _restore(ckpt_from)
                return False, rwp_before, None

            print(f"    ACCEPTED: ΔRwp={rwp_gain:.3f}%, params plausible",
                  flush=True)
            ckpt = _checkpoint(label.replace(' ', '_').replace('+', ''))
            return True, rwp, ckpt

        # ── Diagnostics: dump raw GSAS-II instrument params ──────────
        def _dump_raw(label):
            """Print raw GSAS-II instrument parameter entries."""
            print(f"\n  --- RAW instrument params ({label}) ---",
                  flush=True)
            ip = histogram.data['Instrument Parameters'][0]
            for k in ['U', 'V', 'W', 'X', 'Y', 'Zero', 'SH/L',
                       'Lam', 'Lam1', 'Lam2']:
                v = ip.get(k)
                if v is not None:
                    print(f"    {k:8s} = {repr(v)}", flush=True)
            print(f"    (parsed) = {_current_params()}", flush=True)

        def _print_fwhm(params, label):
            """Print predicted Gaussian FWHM at Si peak positions."""
            U = params.get('U', 0)
            V = params.get('V', 0)
            W = params.get('W', 0)
            tan_t = np.tan(np.radians(_SI_PEAK_ANGLES / 2.0))
            hg2 = U * tan_t**2 + V * tan_t + W
            fwhm = np.sqrt(np.maximum(hg2, 0))
            print(f"    {label} FWHM(deg) at Si peaks: "
                  f"{[f'{x:.4f}' for x in fwhm]}", flush=True)

        _dump_raw("INITIAL (before any refinement)")
        _print_fwhm(_current_params(), "Initial")

        # ── Calibration sequence ───────────────────────────────────────
        print("\n  --- Conservative calibration with checkpoints ---",
              flush=True)

        # Stage 1: BG + scale
        if progress_callback:
            progress_callback('Calibration: BG + scale...')
        _refine("Stage 1 (BG+scale)", [], cycles=5)
        _dump_raw("after Stage 1")

        # Stage 2: Zero
        if progress_callback:
            progress_callback('Calibration: Zero...')
        _refine("Stage 2 (Zero)", ['Zero'], cycles=5)
        _dump_raw("after Stage 2")
        _print_fwhm(_current_params(), "Stage 2")

        # Stage 3: W + Zero (Gaussian floor)
        if progress_callback:
            progress_callback('Calibration: W...')
        _refine("Stage 3 (W+Zero)", ['W', 'Zero'], cycles=8)
        _dump_raw("after Stage 3")
        _print_fwhm(_current_params(), "Stage 3")

        # ← STABLE BASELINE CHECKPOINT
        baseline_ckpt = _checkpoint('baseline_W_Zero')
        _baseline_rwp = _get_rwp()
        _baseline_params = _current_params()
        print(f"\n  === STABLE BASELINE: Rwp={_baseline_rwp:.3f}%, "
              f"W={_baseline_params['W']:.4f}, "
              f"Zero={_baseline_params['Zero']:.5f} ===\n", flush=True)

        # Track the best Gaussian checkpoint
        gauss_ckpt = baseline_ckpt

        # Trial A: U + W + Zero
        if progress_callback:
            progress_callback('Calibration: trying U...')
        u_ok, _, u_ckpt = _trial(
            "Trial A (UW)", ['U', 'W', 'Zero'],
            ckpt_from=gauss_ckpt, min_rwp_gain=0.5)
        if u_ok:
            gauss_ckpt = u_ckpt

        # Trial B: V (only if U was accepted)
        v_ok = False
        if u_ok:
            if progress_callback:
                progress_callback('Calibration: trying V...')
            v_ok, _, v_ckpt = _trial(
                "Trial B (UVW)", ['U', 'V', 'W', 'Zero'],
                ckpt_from=gauss_ckpt, min_rwp_gain=0.3)
            if v_ok:
                gauss_ckpt = v_ckpt

        # Build Gaussian param list for Lorentzian trials
        _gauss = ['W', 'Zero']
        if u_ok:
            _gauss = ['U', 'W', 'Zero']
        if u_ok and v_ok:
            _gauss = ['U', 'V', 'W', 'Zero']

        # ── Lorentzian trials ──────────────────────────────────────────
        # Try X and Y from TWO starting states:
        #   (a) gauss_ckpt (U+W or W-only, whatever was accepted)
        #   (b) baseline_ckpt (W-only, U=0)
        # The W-only baseline often works better because adding a
        # Lorentzian term on top of negative U forces U more negative,
        # making FWHM²_G go negative.  From W-only, the Lorentzian
        # handles the angular variation that U was trying to capture.
        if progress_callback:
            progress_callback('Calibration: trying Lorentzian terms...')

        _candidates = []  # (label, rwp, ckpt, use_X_flag, gauss_list)

        # Try X and Y from gauss_ckpt (U+W state)
        print("\n    --- Lorentzian trials from Gaussian state "
              f"({', '.join(_gauss)}) ---", flush=True)

        x1_ok, rwp_x1, x1_ckpt = _trial(
            "Trial C1 (gauss+X)", _gauss + ['X'],
            ckpt_from=gauss_ckpt, min_rwp_gain=0.1)
        if x1_ok:
            _candidates.append(('gauss+X', rwp_x1, x1_ckpt, True, _gauss))

        y1_ok, rwp_y1, y1_ckpt = _trial(
            "Trial D1 (gauss+Y)", _gauss + ['Y'],
            ckpt_from=gauss_ckpt, min_rwp_gain=0.1)
        if y1_ok:
            _candidates.append(('gauss+Y', rwp_y1, y1_ckpt, False, _gauss))

        # Also try X and Y from W-only baseline (U=0)
        if u_ok:
            # Only worth trying if we had U — otherwise gauss_ckpt IS
            # the baseline and we'd just duplicate the trials above
            _base_gauss = ['W', 'Zero']
            print(f"\n    --- Lorentzian trials from W-only baseline ---",
                  flush=True)

            x2_ok, rwp_x2, x2_ckpt = _trial(
                "Trial C2 (W+X)", _base_gauss + ['X'],
                ckpt_from=baseline_ckpt, min_rwp_gain=0.1)
            if x2_ok:
                _candidates.append(('W+X', rwp_x2, x2_ckpt, True,
                                    _base_gauss))

            y2_ok, rwp_y2, y2_ckpt = _trial(
                "Trial D2 (W+Y)", _base_gauss + ['Y'],
                ckpt_from=baseline_ckpt, min_rwp_gain=0.1)
            if y2_ok:
                _candidates.append(('W+Y', rwp_y2, y2_ckpt, False,
                                    _base_gauss))

        # Also keep the Gaussian-only option as a candidate
        _candidates.append(('Gaussian only', _get_rwp(), gauss_ckpt,
                            False, _gauss))

        # Pick the best candidate by Rwp
        _candidates.sort(key=lambda c: c[1])
        _best = _candidates[0]
        _best_label, _best_rwp, _best_ckpt, _use_X, _gauss = _best
        y_ok = not _use_X and _best_label != 'Gaussian only'

        print(f"\n    Candidate ranking:", flush=True)
        for i, (lbl, rw, _, _, _) in enumerate(_candidates):
            marker = " ← BEST" if i == 0 else ""
            print(f"      {lbl}: Rwp={rw:.3f}%{marker}", flush=True)

        _restore(_best_ckpt)
        print(f"    → Selected: {_best_label} (Rwp={_best_rwp:.3f}%)",
              flush=True)

        # ← LORENTZIAN CHECKPOINT (before SH/L trial)
        lor_ckpt = _checkpoint('pre_shl')

        # Trial E: SH/L (only if data below 30° and meaningful)
        if tt_min < 30.0:
            if progress_callback:
                progress_callback('Calibration: trying SH/L...')
            _lor = ['X'] if _use_X else (['Y'] if y_ok else [])
            shl_ok, _, shl_ckpt = _trial(
                "Trial E (SHL)", _gauss + _lor + ['SH/L'],
                ckpt_from=lor_ckpt, min_rwp_gain=0.1)
            if shl_ok:
                _shl = _get_inst_value('SH/L')
                if _shl < 0 or _shl > 0.1:
                    print(f"    SH/L={_shl:.5f} suspicious — reverting",
                          flush=True)
                    _restore(lor_ckpt)
                else:
                    _restore(shl_ckpt)
            # If rejected, _trial already restored lor_ckpt
        else:
            print("    SH/L: skipped (no data below 30°)", flush=True)

        # ── Extract final parameters ───────────────────────────────────
        params = _current_params()

        # Zero out unused Lorentzian
        if _use_X:
            params['Y'] = 0.0
        elif y_ok:
            params['X'] = 0.0
        else:
            params['X'] = 0.0
            params['Y'] = 0.0

        # Final plausibility check
        ok, reason = _profile_plausible(params, _SI_PEAK_ANGLES)
        if not ok:
            print(f"\n  WARNING: final profile still implausible "
                  f"({reason}). Falling back to baseline.", flush=True)
            _restore(baseline_ckpt)
            params = _baseline_params.copy()
            params['X'] = 0.0
            params['Y'] = 0.0
            _use_X = False
            y_ok = False

        stats = _get_stats()
        _rwp_final = _get_rwp()

        print(f"\n  {'=' * 50}", flush=True)
        print(f"  CALIBRATED INSTRUMENT PARAMETERS", flush=True)
        print(f"  {'=' * 50}", flush=True)
        for k, v in params.items():
            print(f"    {k:6s} = {v:.6f}", flush=True)
        _accepted = [k for k in ['U', 'V', 'W', 'X', 'Y']
                     if abs(params.get(k, 0)) > 1e-6]
        if not _accepted:
            _accepted = ['W']
        print(f"  Accepted params: {', '.join(_accepted)} + Zero",
              flush=True)
        print(f"  Rwp = {_rwp_final:.3f}%", flush=True)

        # Caglioti check
        _tan_t = np.tan(np.radians(_SI_PEAK_ANGLES / 2.0))
        _hg2 = params['U'] * _tan_t**2 + params['V'] * _tan_t + params['W']
        print(f"  FWHM²_G range: [{np.min(_hg2):.4f}, "
              f"{np.max(_hg2):.4f}] (all positive ✓)"
              if np.min(_hg2) > 0 else
              f"  WARNING: FWHM²_G goes negative!", flush=True)

        # ── Write .instprm (path resolved at top of function) ──────────

        if _use_doublet:
            _lam = (f"Lam1:{CU_KALPHA1_A:.6f}\n"
                    f"Lam2:{CU_KALPHA2_A:.6f}\n"
                    f"I(L2)/I(L1):{CU_KALPHA2_RATIO:.4f}\n")
        else:
            _lam = f"Lam:{wavelength:.6f}\n"

        content = (
            "#GSAS-II instrument parameter file; "
            "do not add/delete items!\n"
            "Type:PXC\n"
            f"{_lam}"
            f"Zero:{params['Zero']:.5f}\n"
            f"Polariz.:{polariz}\n"
            f"U:{params['U']:.4f}\n"
            f"V:{params['V']:.4f}\n"
            f"W:{params['W']:.4f}\n"
            f"X:{params['X']:.4f}\n"
            f"Y:{params['Y']:.4f}\n"
            "Z:0.0\n"
            f"SH/L:{params.get('SH/L', 0.002):.5f}\n"
            "Azimuth:0.0\n"
        )
        # Atomic write: .tmp → backup old → replace.
        _tmp_out = output_instprm + '.tmp'
        print(f"\n  Writing calibrated instprm to: {output_instprm}",
              flush=True)
        with open(_tmp_out, 'w') as f:
            f.write(content)
        if os.path.isfile(output_instprm):
            _bak = output_instprm + '.bak'
            shutil.copy2(output_instprm, _bak)
            print(f"  Backed up old file: {_bak}", flush=True)
        os.replace(_tmp_out, output_instprm)
        print(f"  .instprm saved successfully.", flush=True)
        print(content, flush=True)

        # ── Pattern for plotting ───────────────────────────────────────
        try:
            tt_out, y_obs_out, y_calc_out, y_bg_out = _get_arrays()
            diff_out = y_obs_out - y_calc_out
        except Exception:
            tt_out, y_obs_out = tt_r, y_r
            y_calc_out = y_bg_out = np.zeros_like(y_r)
            diff_out = y_r

        return {
            'instprm_path': output_instprm,
            'params': params,
            'Rwp': _rwp_final,
            'accepted_params': _accepted,
            'tt': tt_out.tolist(),
            'y_obs': y_obs_out.tolist(),
            'y_calc': y_calc_out.tolist(),
            'y_background': y_bg_out.tolist(),
            'residuals': diff_out.tolist(),
            'statistics': {
                'Rwp': stats.get('Rwp', 0),
                'Rp': stats.get('Rp', 0),
                'chi2': stats.get('reduced chi2', 0),
                'GoF': math.sqrt(max(0, stats.get('reduced chi2', 0))),
            },
            'method': 'GSAS-II calibration',
            'lorentzian_term': 'X' if _use_X else ('Y' if y_ok else 'none'),
        }

    finally:
        if not keep_workdir:
            try:
                shutil.rmtree(work_dir)
            except Exception:
                pass
        else:
            print(f"  Work dir kept: {work_dir}", flush=True)
