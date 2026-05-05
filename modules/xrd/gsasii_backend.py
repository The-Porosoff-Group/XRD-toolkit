"""
modules/xrd/gsasii_backend.py
GSAS-II integration for Rietveld/Le Bail refinement via GSASIIscriptable.

Requires GSAS-II installed in the Python environment.
Install via:  git clone https://github.com/AdvancedPhotonSource/GSAS-II.git
              cd GSAS-II && pip install .
This module wraps GSASIIscriptable to provide a refinement backend compatible
with the toolkit's result format (same keys as run_lebail / run_rietveld).

SEPARATION RULE (refinement vs. display):
  This module owns the CIF that GSAS-II refines.  The CIF GSAS-II reads
  must be internally consistent: cell parameters and atom positions
  must travel together in the same crystallographic setting.  CifWriter
  may transform an asymmetric unit into a different setting (e.g. axis
  permutation for orthorhombic phases); when that happens the entire
  CIF block — cell + sites + symmetry ops — must be used together, NOT
  combined with cell parameters from the search-row metadata.

  The preview/tick generator (cod_api.get_stick_pattern) is allowed to
  produce a different reflection list for display.  That display
  decision MUST NOT feed back into the CIF construction here.

  Validation gate (post-add_phase, before any refinement) compares
  expected SG/cell/atom-count against what GSAS-II actually imported
  and warns on mismatch.  This catches the W2C-as-P1 class of bug.

  Architecture spec: see CLAUDE OUTPUTS /
  Catalysis-Toolkit_Architecture_v1.md.
"""

import math, os, re, sys, tempfile, warnings
import numpy as np


class _IsolationSkipped(Exception):
    """Sentinel: phase isolation intentionally disabled."""
    pass


# ── GSAS-II availability check ──────────────────────────────────────────────

_GSASII_AVAILABLE = False
_GSASII_IMPORT_ERROR = None

def _add_gsas2pkg_paths():
    """Add gsas2pkg conda install paths to sys.path if not already present.
    gsas2pkg installs to {conda_prefix}/GSAS-II/ rather than site-packages."""
    prefix = sys.prefix  # e.g. C:\catalysis-toolkit\.conda_env
    candidates = [
        os.path.join(prefix, 'GSAS-II', 'GSASII'),   # for import GSASIIscriptable
        os.path.join(prefix, 'GSAS-II'),               # for import GSASII.GSASIIscriptable
        os.path.join(prefix, 'GSAS-II', 'backcompat'), # backcompat shim
    ]
    for p in candidates:
        if os.path.isdir(p) and p not in sys.path:
            sys.path.insert(0, p)

_add_gsas2pkg_paths()

try:
    # New-style import (pip-installed from GitHub, or gsas2pkg via GSAS-II subpackage)
    import GSASII.GSASIIscriptable as G2sc
    _GSASII_AVAILABLE = True
except ImportError:
    try:
        # Direct import (gsas2pkg installs GSASII/ dir; backcompat or GSASII on path)
        import GSASIIscriptable as G2sc
        _GSASII_AVAILABLE = True
    except ImportError as e:
        _GSASII_IMPORT_ERROR = str(e)

try:
    from .crystallography import (
        compute_fit_statistics, cell_volume, molar_mass_from_formula,
        tch_fwhm_eta, size_from_Y, scherrer_size,
        generate_reflections, parse_cif, caglioti_fwhm,
        expand_sites_from_cif,
    )
except ImportError:
    from crystallography import (
        compute_fit_statistics, cell_volume, molar_mass_from_formula,
        tch_fwhm_eta, size_from_Y, scherrer_size,
        generate_reflections, parse_cif, caglioti_fwhm,
        expand_sites_from_cif,
    )


def is_available():
    """Return True if GSASIIscriptable is importable and functional."""
    if not _GSASII_AVAILABLE:
        return False
    # Quick functional check: verify G2Project is callable
    # (catches partial installs where Fortran extensions are missing)
    try:
        if not hasattr(G2sc, 'G2Project'):
            return False
    except Exception:
        return False
    return True


def import_error():
    """Return the import error message, or None if available."""
    return _GSASII_IMPORT_ERROR


# ─────────────────────────────────────────────────────────────────────────────
# SPACE GROUP TABLE  (number → H-M symbol for CIF generation)
# ─────────────────────────────────────────────────────────────────────────────

_SG_HM = {
    # Triclinic / Monoclinic
    1: 'P 1', 2: 'P -1', 3: 'P 2', 4: 'P 21', 5: 'C 2',
    6: 'P m', 7: 'P c', 8: 'C m', 9: 'C c', 10: 'P 2/m',
    11: 'P 21/m', 12: 'C 2/m', 13: 'P 2/c', 14: 'P 21/c', 15: 'C 2/c',
    # Orthorhombic (common catalyst phases)
    16: 'P 2 2 2', 19: 'P 21 21 21', 26: 'P m c 21', 29: 'P c a 21',
    31: 'P m n 21', 33: 'P n a 21', 36: 'C m c 21', 40: 'A m a 2',
    44: 'I m m 2', 46: 'I m a 2',
    47: 'P m m m', 51: 'P m m a', 55: 'P b a m', 57: 'P b c m',
    58: 'P n n m', 59: 'P m m n', 60: 'P b c n', 61: 'P b c a',
    62: 'P n m a', 63: 'C m c m', 64: 'C m c a', 65: 'C m m m',
    66: 'C c c m', 69: 'F m m m', 70: 'F d d d', 71: 'I m m m',
    72: 'I b a m', 74: 'I m m a',
    # Tetragonal
    75: 'P 4', 82: 'I -4', 83: 'P 4/m', 84: 'P 42/m',
    85: 'P 4/n', 87: 'I 4/m', 88: 'I 41/a',
    99: 'P 4 m m', 107: 'I 4 m m', 115: 'P -4 m 2',
    119: 'I -4 m 2', 121: 'I -4 2 m', 122: 'I -4 2 d',
    123: 'P 4/m m m', 129: 'P 4/n m m', 131: 'P 42/m m c',
    136: 'P 42/m n m', 139: 'I 4/m m m', 140: 'I 4/m c m',
    141: 'I 41/a m d', 142: 'I 41/a c d',
    # Trigonal / Hexagonal
    143: 'P 3', 146: 'R 3', 147: 'P -3', 148: 'R -3',
    150: 'P 3 2 1', 152: 'P 31 2 1', 154: 'P 32 2 1',
    155: 'R 3 2', 156: 'P 3 m 1', 157: 'P 3 1 m',
    160: 'R 3 m', 161: 'R 3 c', 162: 'P -3 1 m', 163: 'P -3 1 c',
    164: 'P -3 m 1', 165: 'P -3 c 1', 166: 'R -3 m', 167: 'R -3 c',
    168: 'P 6', 173: 'P 63', 174: 'P -6', 175: 'P 6/m', 176: 'P 63/m',
    183: 'P 6 m m', 186: 'P 63 m c', 187: 'P -6 m 2', 189: 'P -6 2 m',
    191: 'P 6/m m m', 193: 'P 63/m c m', 194: 'P 63/m m c',
    # Cubic
    195: 'P 2 3', 196: 'F 2 3', 197: 'I 2 3', 198: 'P 21 3',
    199: 'I 21 3', 200: 'P m -3', 201: 'P n -3', 202: 'F m -3',
    203: 'F d -3', 204: 'I m -3', 205: 'P a -3', 206: 'I a -3',
    207: 'P 4 3 2', 209: 'F 4 3 2', 211: 'I 4 3 2',
    212: 'P 43 3 2', 213: 'P 41 3 2', 214: 'I 41 3 2',
    215: 'P -4 3 m', 216: 'F -4 3 m', 217: 'I -4 3 m',
    218: 'P -4 3 n', 219: 'F -4 3 c', 220: 'I -4 3 d',
    221: 'P m -3 m', 222: 'P n -3 n', 223: 'P m -3 n',
    224: 'P n -3 m', 225: 'F m -3 m', 226: 'F m -3 c',
    227: 'F d -3 m', 228: 'F d -3 c', 229: 'I m -3 m',
    230: 'I a -3 d',
}


def _get_hm_symbol(sg_num):
    """Get H-M symbol for a space group number.

    Uses the static table first, then tries pymatgen as a fallback.
    """
    if sg_num in _SG_HM:
        return _SG_HM[sg_num]
    try:
        from pymatgen.symmetry.groups import SpaceGroup
        sg = SpaceGroup.from_int_number(sg_num)
        return sg.symbol
    except Exception:
        pass
    warnings.warn(f"Space group {sg_num} not in lookup table — using 'P 1'. "
                  f"Add it to _SG_HM in gsasii_backend.py for correct results.")
    return 'P 1'


def _hm_symbol_to_number(symbol):
    """Best-effort H-M symbol to International Tables number."""
    text = str(symbol or '').strip().strip("'\"")
    if not text:
        return 0
    try:
        return int(text)
    except (TypeError, ValueError):
        pass

    def _norm(s):
        return re.sub(r'[\s_]+', '', str(s)).lower()

    target = _norm(text)
    for num, hm in _SG_HM.items():
        if _norm(hm) == target:
            return int(num)
    try:
        from pymatgen.symmetry.groups import SpaceGroup
        return int(SpaceGroup(text).int_number)
    except Exception:
        return 0


# ── Instrument profiles ─────────────────────────────────────────────────
# Named profiles bundle geometry, displacement model, polarization, SH/L,
# preferred-orientation default, and zero seed for each known instrument.
# These are provisional — a measured .instprm file ALWAYS takes priority
# over profile defaults for U/V/W/X/Y/SH/L.
#
# Usage:
#   1. Auto-detect from filename/header/metadata via infer_instrument().
#   2. Override explicitly via params['instrument'] = 'smartlab'.
#   3. Fall back to DEFAULT_INSTRUMENT if inference fails.
INSTRUMENT_PROFILES = {
    'synergy_s': {
        'label':    'Synergy-S capillary/transmission',
        'geometry': 'capillary',
        'displacement_param': 'DisplaceY',
        'zero_seed': -0.25,    # from Synergy-Dualflex .par (WC/W2C config)
        'polariz':  0.5,
        'sh_l':     0.002,
        'preferred_orientation_default': 'off',
        'sigma_inflation_K': 5.0,  # 2D-integrated σ underestimates uncertainty
        'instprm_filename': 'synergy_s_Si640g.instprm',
        'notes': 'Measured from NIST SRM 640g Si standard.',
    },
    'smartlab': {
        'label':    'Rigaku SmartLab flat plate / Bragg-Brentano',
        'geometry': 'bragg_brentano',
        'displacement_param': 'DisplaceX',
        'zero_seed': -0.027,   # measured from NIST Si 640g
        'polariz':  0.7,       # SmartLab manual calibration profile
        'sh_l':     0.002,     # SmartLab manual calibration profile
        'calibration_allow_x': False,
        'calibration_allow_y': False,
        'calibration_refine_sh_l': False,
        'calibration_fixed_sh_l': 0.002,
        'calibration_u_min_rwp_gain': 0.0,
        'calibration_v_min_rwp_gain': 0.0,
        'preferred_orientation_default': 'auto',
        'sigma_inflation_K': 1.0,  # σ ≈ √I is already correct for BB
        'instprm_filename': 'smartlab_Si640g.instprm',  # auto-locate
        'notes': 'Measured from NIST SRM 640g Si standard, 2026-04-30.',
    },
}
DEFAULT_INSTRUMENT = 'smartlab'

# Legacy defaults — kept for backward compat and as fallback values.
# When an instrument profile is active, these are overridden by profile
# values.
DEFAULT_POLARIZ = 0.5
DEFAULT_SH_L = 0.002
DEFAULT_U = 2.0            # Caglioti U initial guess (centideg²) — refined
DEFAULT_V = -2.0           # Caglioti V initial guess (centideg²) — refined
DEFAULT_W = 5.0            # Caglioti W initial guess (centideg²) — refined
DEFAULT_X = 0.0            # Lorentzian X initial guess (centideg) — refined
DEFAULT_Y = 0.0            # Lorentzian Y initial guess (centideg) — refined
DEFAULT_B_ISO = 0.5        # Fallback B_iso (Å²) when GSAS-II extraction fails


def infer_instrument(filepath=None, metadata=None, raw_header=None):
    """Infer instrument identity from filename, metadata, or header text.

    Returns (instrument_key, reason_string).  Falls back to
    DEFAULT_INSTRUMENT if nothing matches.
    """
    metadata = metadata or {}
    text = ' '.join([
        str(filepath or ''),
        str(metadata),
        str(raw_header or ''),
    ]).lower()

    # Strong filename/header hints
    if 'synergy' in text or 'crysalis' in text or 'dualflex' in text:
        return 'synergy_s', 'filename/header contains Synergy/CrysAlis/Dualflex'
    if 'smartlab' in text or 'rigaku' in text:
        return 'smartlab', 'filename/header contains SmartLab/Rigaku'

    # Weak format hints
    fmt = str(metadata.get('format', '')).lower()
    if fmt == 'stepscan':
        return 'smartlab', 'weak inference from StepScan format'
    if fmt == 'powdergraph':
        return 'synergy_s', 'weak inference from PowderGraph/integrated format'

    return DEFAULT_INSTRUMENT, (
        f'default fallback ({DEFAULT_INSTRUMENT}); '
        f'instrument not identifiable from file')


def _get_expanded_sites(cif_text, spacegroup_number=None):
    """Get symmetry-expanded unit-cell sites for structure factor computation.

    Uses the shared expand_sites_from_cif (pymatgen) first, then falls back
    to raw parse_cif (may be asymmetric unit only).
    """
    if not cif_text:
        return None
    sites = expand_sites_from_cif(cif_text)
    if sites:
        return sites
    # Fallback: raw parse_cif (may be asymmetric unit only)
    try:
        parsed = parse_cif(cif_text)
        return parsed.get('sites') or None
    except Exception:
        return None


def _cif_already_has_asymmetric_unit(cif_text, declared_sg=None):
    """Check if a CIF file already contains asymmetric-unit sites (not full cell).

    CIF files from COD and from pymatgen's CifWriter(symprec=...) typically
    list only the asymmetric unit plus symmetry operations.  In contrast,
    CIF files generated by dumping a pymatgen Structure without symmetry
    detection list the full unit cell with space group P1.

    Detection strategy:
      1. Parse sites from the raw CIF text (what is literally written)
      2. Require matching SG plus symmetry-operation tags
      3. Treat that atom loop as the asymmetric unit without asking
         pymatgen to re-expand it; CifParser can exit the interpreter for
         some otherwise valid CIFs.

    Returns (True, raw_sites) if already asymmetric, (False, None) otherwise.
    """
    if not cif_text or not declared_sg or declared_sg <= 1:
        return False, None

    try:
        # Count sites literally listed in the CIF
        raw_parsed = parse_cif(cif_text)
        raw_sites = raw_parsed.get('sites') or []
        raw_sg = raw_parsed.get('spacegroup_number', 1)
        if not raw_sites or raw_sg <= 1:
            return False, None

        # Check if the CIF contains symmetry operations (loop_ _symmetry_equiv
        # or _space_group_symop) — this is the hallmark of an asymmetric-unit CIF
        has_symops = ('_symmetry_equiv_pos_as_xyz' in cif_text or
                      '_space_group_symop_operation_xyz' in cif_text or
                      '_space_group_symop.operation_xyz' in cif_text)

        if not has_symops:
            return False, None

        # Check that the CIF's own SG is compatible with the declared SG.
        # If CifWriter detected a DIFFERENT SG (e.g. SG 25 instead of 60),
        # the sites are reduced for the wrong symmetry and must NOT be used
        # with the declared SG — that would cause GSAS-II to expand with
        # wrong symmetry operations.
        if raw_sg != declared_sg:
            print(f"  CIF declares SG {raw_sg} but phase expects SG "
                  f"{declared_sg} — cannot use CIF sites as-is "
                  f"(SG mismatch)", flush=True)
            return False, None

        print(f"  CIF already contains symmetry operations: "
              f"{len(raw_sites)} atom site(s), SG {raw_sg}. Treating "
              f"the atom loop as the asymmetric unit.", flush=True)
        return True, raw_sites

    except BaseException as e:
        print(f"  Warning: asymmetric-unit detection failed: {e}", flush=True)

    return False, None


def _cifwriter_asymmetric_unit(cif_text, declared_sg=None,
                               return_full_cif=False):
    """Use pymatgen CifWriter to generate a proper asymmetric-unit CIF.

    CifWriter with symprec is purpose-built for reducing a full-cell
    structure to the asymmetric unit with correct symmetry operations.
    This is more robust than manually using SpacegroupAnalyzer because
    CifWriter handles the full reduction pipeline internally, including
    site merging and Wyckoff position assignment.

    IMPORTANT: CifWriter may transform the structure to the standard ITA
    setting, which can permute axes relative to the original structure.
    When return_full_cif=True, the complete CIF text is returned so the
    caller gets atom positions AND cell parameters in the SAME convention.
    Mixing CifWriter positions with original cell params causes axis
    mismatch → wrong structure factors → zero intensity in Rietveld.

    Returns:
        If return_full_cif=False:
            list of (element, x, y, z, occupancy) tuples, or None
        If return_full_cif=True:
            (sites_list, full_cif_text) tuple, or (None, None)
    """
    if not cif_text:
        return None

    try:
        from pymatgen.io.cif import CifParser, CifWriter
        from pymatgen.core import Structure as PmgStruct

        # Parse the input CIF to get the full conventional cell
        try:
            parser = CifParser.from_str(cif_text, occupancy_tolerance=100.0)
        except (AttributeError, TypeError):
            import tempfile as _tf, os as _os
            _fd, _tmp = _tf.mkstemp(suffix='.cif')
            with _os.fdopen(_fd, 'w') as _f:
                _f.write(cif_text)
            parser = CifParser(_tmp, occupancy_tolerance=100.0)
            _os.unlink(_tmp)
        structs = parser.parse_structures(primitive=False)
        if not structs:
            return None

        struct0 = structs[0]
        full_cell_n = len(struct0)

        # Try CifWriter with multiple symprec values.
        # CifWriter internally uses SpacegroupAnalyzer but goes through
        # the full CIF generation pipeline which handles edge cases better
        # than raw SpacegroupAnalyzer.get_symmetrized_structure().
        for symprec in (0.01, 0.05, 0.1, 0.2):
            try:
                writer = CifWriter(struct0, symprec=symprec)
                # Extract the reduced CIF text
                import tempfile as _tf2, os as _os2
                _fd2, _tmp2 = _tf2.mkstemp(suffix='.cif')
                with _os2.fdopen(_fd2, 'w') as _f2:
                    pass  # just create the file
                writer.write_file(_tmp2)
                with open(_tmp2, 'r') as _f2:
                    reduced_cif = _f2.read()
                _os2.unlink(_tmp2)

                # Parse the reduced CIF to get the asymmetric unit sites
                reduced_parsed = parse_cif(reduced_cif)
                reduced_sites = reduced_parsed.get('sites') or []
                reduced_sg = reduced_parsed.get('spacegroup_number', 1)

                if not reduced_sites:
                    continue

                # Validate: reduced must be smaller than full cell
                if declared_sg and declared_sg > 1:
                    if len(reduced_sites) >= full_cell_n:
                        print(f"  CifWriter(symprec={symprec}): "
                              f"{len(reduced_sites)} sites = full cell "
                              f"({full_cell_n}) — no reduction, skipping",
                              flush=True)
                        continue

                # Validate: re-expand the reduced CIF and check site count
                try:
                    parser2 = CifParser.from_str(
                        reduced_cif, occupancy_tolerance=100.0)
                    structs2 = parser2.parse_structures(primitive=False)
                    if structs2:
                        reexpanded_n = len(structs2[0])
                        if abs(reexpanded_n - full_cell_n) > 1:
                            print(f"  CifWriter(symprec={symprec}): "
                                  f"re-expansion gave {reexpanded_n} sites, "
                                  f"expected {full_cell_n} — skipping",
                                  flush=True)
                            continue
                except Exception:
                    pass  # skip validation, still try the sites

                # Validate: detected SG should be compatible with declared
                if (declared_sg and declared_sg > 1
                        and reduced_sg > 1 and reduced_sg != declared_sg):
                    # Allow if they share the same crystal system at least
                    # (CifWriter might pick a different but valid setting)
                    print(f"  CifWriter(symprec={symprec}): detected "
                          f"SG {reduced_sg} vs declared {declared_sg} "
                          f"— using CifWriter result anyway", flush=True)

                print(f"  CifWriter asymmetric unit: {full_cell_n} → "
                      f"{len(reduced_sites)} sites "
                      f"(symprec={symprec}, SG {reduced_sg})", flush=True)
                if return_full_cif:
                    return reduced_sites, reduced_cif
                return reduced_sites

            except Exception as exc:
                print(f"  CifWriter(symprec={symprec}): failed ({exc})",
                      flush=True)
                continue

    except BaseException as e:
        print(f"  Warning: CifWriter fallback failed: {e}", flush=True)

    if return_full_cif:
        return None, None
    return None


def _reduce_to_asymmetric_unit(cif_text, declared_sg=None):
    """Reduce a full-unit-cell CIF to asymmetric unit for GSAS-II.

    GSAS-II applies symmetry operations itself, so it expects only the
    asymmetric unit in the CIF.  If we give it the full unit cell plus
    the space group, it over-expands and generates wrong reflections.

    Reduction strategy (in order of preference):
      0. Check if CIF already contains asymmetric unit (COD CIFs, or
         Materials Project CIFs written by CifWriter with symprec).
         If so, use those sites directly — no reduction needed.
      1. SpacegroupAnalyzer symmetrized structure (original method,
         works well for simple metals).
      2. CifWriter-based reduction (more robust for complex structures
         like W2C Pbcn where SpacegroupAnalyzer struggles).
      3. Clustering-based reduction using declared space group.
      4. Raw CIF sites fallback (last resort, may cause over-expansion).

    declared_sg: expected space group number from phase dict (used for
                 validation, not to override pymatgen's detection).
    """
    if not cif_text:
        return []

    # ── Step 0: Check if CIF already has asymmetric unit ────────────────
    # CIFs from COD and from CifWriter(symprec=...) already contain the
    # asymmetric unit plus symmetry operations.  Trying to "reduce" these
    # is a no-op at best and destructive at worst (SpacegroupAnalyzer may
    # misidentify equivalent sites in compact structures).
    already_asym, asym_sites = _cif_already_has_asymmetric_unit(
        cif_text, declared_sg)
    if already_asym and asym_sites:
        return asym_sites

    try:
        from pymatgen.io.cif import CifParser
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        try:
            parser = CifParser.from_str(cif_text, occupancy_tolerance=100.0)
        except (AttributeError, TypeError):
            import tempfile as _tf, os as _os
            _fd, _tmp = _tf.mkstemp(suffix='.cif')
            with _os.fdopen(_fd, 'w') as _f:
                _f.write(cif_text)
            parser = CifParser(_tmp, occupancy_tolerance=100.0)
            _os.unlink(_tmp)
        structs = parser.parse_structures(primitive=False)
        if structs:
            full_cell_n = len(structs[0])  # total sites in full unit cell

            # ── Step 1: SpacegroupAnalyzer (original method) ────────────
            # Try tight tolerance first (0.01 Å), fall back to looser (0.1 Å).
            # Tight tolerance avoids merging non-equivalent sites in compact
            # cells like W2C Pbcn where heavy atoms are close together.
            sites = None
            for symprec in (0.01, 0.05, 0.1):
                try:
                    sga = SpacegroupAnalyzer(structs[0], symprec=symprec)
                    detected_sg = sga.get_space_group_number()
                    sym_struct = sga.get_symmetrized_structure()
                    candidate = []
                    for equiv_sites in sym_struct.equivalent_sites:
                        site = equiv_sites[0]
                        frac = site.frac_coords % 1.0
                        el = str(site.specie)
                        occ = 1.0
                        if hasattr(site, 'species'):
                            sp = site.species
                            if hasattr(sp, 'num_atoms'):
                                occ = float(sp.num_atoms)
                        candidate.append((el, float(frac[0]), float(frac[1]),
                                          float(frac[2]), occ))

                    # Validate: the number of equivalent sites summed across
                    # all groups should equal the full unit cell site count.
                    expanded_n = sum(
                        len(eq) for eq in sym_struct.equivalent_sites)

                    # CRITICAL: for SG > 1, the asymmetric unit must be
                    # SMALLER than the full cell. If SpacegroupAnalyzer falls
                    # back to P1, every atom is "unique" and len(candidate)
                    # == full_cell_n — which gives GSAS-II the full cell
                    # labeled as asymmetric unit, causing over-expansion.
                    if declared_sg and declared_sg > 1:
                        if len(candidate) >= full_cell_n:
                            print(f"  symprec={symprec}: detected SG {detected_sg}, "
                                  f"{len(candidate)} unique sites = full cell "
                                  f"({full_cell_n}) — reduction failed, "
                                  f"trying looser tolerance", flush=True)
                            continue
                        # Also check detected SG is compatible (same or higher
                        # symmetry than declared)
                        if detected_sg < declared_sg:
                            print(f"  symprec={symprec}: detected SG {detected_sg} "
                                  f"< declared SG {declared_sg} — trying "
                                  f"looser tolerance", flush=True)
                            continue

                    if expanded_n == full_cell_n and candidate:
                        # ── CRITICAL: Validate representative positions ──
                        # SpacegroupAnalyzer picks equiv_sites[0] as the
                        # representative.  This position may be in a non-
                        # standard setting (e.g. C in W2C Pbcn placed at
                        # (0.5,0.75,0.878) instead of standard 4c (0,y,¼)).
                        # When GSAS-II applies standard operators to a non-
                        # standard representative, it generates the wrong
                        # number of atoms (e.g. 8C instead of 4C → 16 total
                        # instead of 12 → completely wrong structure factors).
                        #
                        # Fix: re-expand with from_spacegroup (which uses
                        # standard ITA operators) and verify atom count.
                        try:
                            from pymatgen.core import Structure as _PmgStruct
                            _test = _PmgStruct.from_spacegroup(
                                detected_sg, structs[0].lattice,
                                [s[0] for s in candidate],
                                [[s[1], s[2], s[3]] for s in candidate],
                                tol=0.01)
                            if abs(len(_test) - full_cell_n) > 1:
                                print(f"  symprec={symprec}: from_spacegroup "
                                      f"gives {len(_test)} atoms (expected "
                                      f"{full_cell_n}) — representatives in "
                                      f"non-standard setting, trying CifWriter "
                                      f"fallback", flush=True)
                                continue
                        except Exception as _fsg_e:
                            print(f"  symprec={symprec}: from_spacegroup "
                                  f"validation failed ({_fsg_e}), continuing "
                                  f"anyway", flush=True)

                        print(f"  Asymmetric unit reduction: {full_cell_n} → "
                              f"{len(candidate)} sites (symprec={symprec}, "
                              f"detected SG {detected_sg})", flush=True)
                        sites = candidate
                        break
                except Exception as exc:
                    print(f"  symprec={symprec}: failed ({exc})", flush=True)
                    continue

            if sites:
                return sites

            # ── Step 2: CifWriter-based reduction ───────────────────────
            # CifWriter handles the full asymmetric-unit extraction pipeline
            # internally and is more robust than raw SpacegroupAnalyzer for
            # complex structures (carbides, oxides with close atom spacing).
            #
            # CRITICAL: CifWriter may transform to the standard ITA setting,
            # which can PERMUTE AXES relative to the original structure.
            # We return the full CIF text so _build_conventional_cif uses
            # matching cell params + atom positions (same axis convention).
            # Without this, cell params from the phase dict (original axes)
            # + atom positions from CifWriter (standard axes) = physical
            # mismatch → wrong structure factors → zero Rietveld intensity.
            print("  SpacegroupAnalyzer reduction failed for all symprec "
                  "values — trying CifWriter fallback...", flush=True)
            cw_sites, cw_full_cif = _cifwriter_asymmetric_unit(
                cif_text, declared_sg, return_full_cif=True)
            if cw_sites and cw_full_cif:
                # Return the full CIF string — _build_conventional_cif
                # will detect the string type and use it directly.
                return cw_full_cif
            if cw_sites:
                return cw_sites

            # ── Step 3: Clustering-based reduction ──────────────────────
            # All symprec values failed — try using declared SG directly
            # with pymatgen's from_spacegroup to rebuild the asymmetric unit
            if declared_sg and declared_sg > 1:
                try:
                    from pymatgen.core import Lattice, Structure as PmgStruct
                    struct0 = structs[0]
                    lat = struct0.lattice
                    # Get unique elements and approximate positions by
                    # clustering sites within the full cell
                    unique_sites = _cluster_to_asymmetric(struct0, declared_sg)
                    if unique_sites and len(unique_sites) < full_cell_n:
                        # Validate by re-expanding
                        test_struct = PmgStruct.from_spacegroup(
                            declared_sg, lat,
                            [s[0] for s in unique_sites],
                            [[s[1], s[2], s[3]] for s in unique_sites])
                        if abs(len(test_struct) - full_cell_n) <= 1:
                            print(f"  Asymmetric unit via clustering: "
                                  f"{full_cell_n} → {len(unique_sites)} sites "
                                  f"(SG {declared_sg})", flush=True)
                            return unique_sites
                except Exception as exc:
                    print(f"  Clustering fallback failed: {exc}", flush=True)

    except BaseException as e:
        print(f"  Warning: pymatgen asymmetric-unit reduction failed: {e}",
              flush=True)

    # ── Step 4: Raw CIF sites fallback ──────────────────────────────────
    # Usually fine for COD CIFs which list asymmetric unit.
    # For MP/pymatgen CIFs which may list the full unit cell, this is
    # DANGEROUS — GSAS-II will over-expand.  We detect this case and
    # fall back to P1 (no symmetry expansion) rather than risk wrong
    # reflections from double-expansion.
    try:
        parsed = parse_cif(cif_text)
        raw_sites = parsed.get('sites') or []
        raw_sg = parsed.get('spacegroup_number', 1)

        if raw_sites and declared_sg and declared_sg > 1:
            # Check if these are likely full-cell sites by comparing
            # site count with what we'd expect for the asymmetric unit.
            # If the CIF has no symmetry operations but declares SG > 1,
            # the sites are probably the full cell (from pymatgen Structure
            # dump without symmetry detection).
            has_symops = ('_symmetry_equiv_pos_as_xyz' in cif_text or
                          '_space_group_symop_operation_xyz' in cif_text or
                          '_space_group_symop.operation_xyz' in cif_text)
            if not has_symops:
                print(f"  WARNING: raw CIF has {len(raw_sites)} sites, "
                      f"SG {declared_sg}, but NO symmetry operations — "
                      f"sites are likely full cell. Returning raw sites "
                      f"with P1 flag to prevent over-expansion.", flush=True)
                # Signal to _build_conventional_cif that these are full-cell
                # sites and should NOT be expanded by the declared SG.
                # We do this by returning them tagged — the caller
                # (_build_conventional_cif) will detect the mismatch and
                # write the CIF with P1 instead.
                return raw_sites
            else:
                print(f"  Using raw CIF sites ({len(raw_sites)} atoms) "
                      f"with SG {declared_sg} — CIF has symmetry operations, "
                      f"sites should be asymmetric unit.", flush=True)
        return raw_sites
    except Exception:
        return []


def _cluster_to_asymmetric(struct, sg_num):
    """Cluster full-cell sites into asymmetric unit using the declared SG.

    For each unique element, tries positions from the full cell one at a
    time, expanding each with the declared space group. Keeps the position
    that best reproduces the observed atom count for that element.

    Returns list of (element, x, y, z, occupancy) or None on failure.
    """
    from pymatgen.core import Lattice, Structure as PmgStruct
    import numpy as _np

    lat = struct.lattice
    # Group sites by element
    element_sites = {}
    for site in struct:
        el = str(site.specie)
        frac = tuple(float(c) % 1.0 for c in site.frac_coords)
        element_sites.setdefault(el, []).append(frac)

    result = []
    used_positions = []  # track what's been accounted for

    for el, positions in element_sites.items():
        target_count = len(positions)
        remaining = list(positions)
        el_asym = []

        while remaining:
            # Try each remaining position as a candidate asymmetric site
            best_pos = None
            best_count = 0
            for pos in remaining:
                try:
                    test = PmgStruct.from_spacegroup(
                        sg_num, lat, [el], [list(pos)])
                    count = sum(1 for s in test if str(s.specie) == el)
                    if count > best_count and count <= len(remaining):
                        best_count = count
                        best_pos = pos
                except Exception:
                    continue

            if best_pos is None or best_count == 0:
                break

            el_asym.append(best_pos)
            # Remove the positions accounted for by this Wyckoff site
            # (find closest matches to the expanded positions)
            try:
                expanded = PmgStruct.from_spacegroup(
                    sg_num, lat, [el], [list(best_pos)])
                for exp_site in expanded:
                    exp_frac = tuple(float(c) % 1.0 for c in exp_site.frac_coords)
                    # Remove the closest match from remaining
                    closest_idx = None
                    closest_dist = 999
                    for i, rem in enumerate(remaining):
                        dist = sum((a - b)**2 for a, b in zip(exp_frac, rem))**0.5
                        # Also check wrapped distances
                        dist_wrap = sum(min((a-b)%1, (b-a)%1)**2
                                       for a, b in zip(exp_frac, rem))**0.5
                        d = min(dist, dist_wrap)
                        if d < closest_dist:
                            closest_dist = d
                            closest_idx = i
                    if closest_idx is not None and closest_dist < 0.1:
                        remaining.pop(closest_idx)
            except Exception:
                break

        for pos in el_asym:
            result.append((el, pos[0], pos[1], pos[2], 1.0))

    return result if result else None


def _cif_policy(ph):
    """Determine CIF handling policy for a phase.

    Returns 'default' for all phases.  The legacy 'mp_w2c_pbcn_compat'
    branch was retired 2026-05-01 after the Step-0 CIF audit
    (see Catalysis-Toolkit_CIF-Audit_v1.md).  W2C now ships as the
    canonical Pbcn fixture (mp_api._LOCAL_FIXTURES['mp-2034']) and
    flows through _build_conventional_cif unchanged.

    The function is kept as a stub so future per-phase routing can be
    added here without re-introducing a fork upstream.
    """
    return 'default'


def _build_conventional_cif(ph):
    """
    Build a synthetic CIF string using the phase dict's (conventional) cell
    parameters and atom sites.

    This ensures GSAS-II always sees a CIF consistent with the conventional
    cell, even when the original CIF used a primitive setting (common with
    Materials Project data).  The space group is written explicitly so that
    GSAS-II applies the correct cell-parameter constraints.

    If the asymmetric-unit reduction fails for a complex structure, this
    function detects the failure and returns the ORIGINAL CIF text instead
    of building a synthetic CIF with full-cell sites + the declared space
    group (which would cause GSAS-II to double-expand).

    NOTE: Phases that need a specific canonical CIF (e.g. mp-2034 W2C)
    should ship that CIF through the phase dict's 'cif_text' field — see
    mp_api._LOCAL_FIXTURES.  The old mp_w2c_pbcn_compat early-return was
    retired 2026-05-01 after the Step-0 audit.
    """
    a   = ph.get('a', 4.0)
    b   = ph.get('b', a)
    c   = ph.get('c', a)
    al  = ph.get('alpha', 90.0)
    be  = ph.get('beta', 90.0)
    ga  = ph.get('gamma', 90.0)
    sg  = ph.get('spacegroup_number', 1)
    formula = ph.get('formula', '')
    Z   = ph.get('Z', '')

    # H-M symbol — try phase dict first, then dynamic lookup
    hm = ph.get('spacegroup', '') or ph.get('spacegroup_name', '')
    if not hm:
        hm = _get_hm_symbol(sg)

    # Get atom sites reduced to the asymmetric unit.
    # GSAS-II applies space-group symmetry itself, so we must NOT give
    # it the full unit cell — that would cause double-counting and wrong
    # reflections (e.g. ghost peaks at ~25° for A15-W).
    cif_text = ph.get('cif_text', '')
    if cif_text:
        already_asym, _asym_sites = _cif_already_has_asymmetric_unit(
            cif_text, declared_sg=sg)
        if already_asym:
            print(f"  _build_conventional_cif [{formula or '?'}]: using "
                  f"source CIF directly (matching SG/cell/sites).",
                  flush=True)
            return cif_text
    sites = _reduce_to_asymmetric_unit(cif_text, declared_sg=sg) if cif_text else []

    # ── CifWriter full CIF passthrough ────────────────────────────────
    # When CifWriter was used for asymmetric-unit reduction, it may have
    # transformed the structure to the standard ITA setting (permuting
    # axes). In this case _reduce_to_asymmetric_unit returns the FULL
    # CIF text (a string) rather than a sites list.  Using this CIF
    # directly ensures atom positions and cell parameters are in the
    # SAME axis convention.  If we instead extracted just the positions
    # and combined them with the original phase dict cell params, axis
    # mismatch would put atoms at wrong physical locations → zero
    # structure factors → phase scale collapses in Rietveld.
    if isinstance(sites, str):
        print(f"  _build_conventional_cif [{formula or '?'}]: using full "
              f"CifWriter CIF directly (axis-consistent cell + positions)",
              flush=True)
        return sites

    # ── Expansion validation: verify that GSAS-II will generate the ──
    # correct number of atoms from the asymmetric unit + space group.
    # This catches cases where the representative positions are in a
    # non-standard setting (e.g. SpacegroupAnalyzer returning C at
    # (0.5,0.75,0.878) instead of standard 4c (0,y,0.25) for Pbcn).
    if sites and sg > 1 and cif_text:
        try:
            from pymatgen.core import Structure as _PmgStruct
            from pymatgen.core import Lattice as _PmgLat
            _lat = _PmgLat.from_parameters(a, b, c, al, be, ga)
            _test = _PmgStruct.from_spacegroup(
                sg, _lat,
                [s[0] for s in sites],
                [[s[1], s[2], s[3]] for s in sites],
                tol=0.01)
            # Get expected full cell count from the original CIF
            from pymatgen.io.cif import CifParser as _CP
            try:
                _p = _CP.from_str(cif_text, occupancy_tolerance=100.0)
            except (AttributeError, TypeError):
                import tempfile as _tf, os as _os
                _fd, _tmp = _tf.mkstemp(suffix='.cif')
                with _os.fdopen(_fd, 'w') as _f:
                    _f.write(cif_text)
                _p = _CP(_tmp, occupancy_tolerance=100.0)
                _os.unlink(_tmp)
            _structs = _p.parse_structures(primitive=False)
            _expected_n = len(_structs[0]) if _structs else len(_test)

            if abs(len(_test) - _expected_n) > 1:
                print(f"  _build_conventional_cif: from_spacegroup gives "
                      f"{len(_test)} atoms (expected {_expected_n}) — "
                      f"atom positions may be in non-standard setting. "
                      f"Attempting CifWriter rescue...", flush=True)
                # Try CifWriter to get correct standard-setting CIF
                from pymatgen.io.cif import CifWriter as _CW
                _struct0 = _structs[0]
                for _symprec in (0.01, 0.05, 0.1, 0.2):
                    try:
                        _writer = _CW(_struct0, symprec=_symprec)
                        import tempfile as _tf3, os as _os3
                        _fd3, _tmp3 = _tf3.mkstemp(suffix='.cif')
                        with _os3.fdopen(_fd3, 'w') as _f3:
                            pass
                        _writer.write_file(_tmp3)
                        with open(_tmp3, 'r') as _f3:
                            _rescued_cif = _f3.read()
                        _os3.unlink(_tmp3)
                        # Validate the rescued CIF
                        _rp = parse_cif(_rescued_cif)
                        _rsites = _rp.get('sites') or []
                        _rsg = _rp.get('spacegroup_number', 1)
                        if _rsites and len(_rsites) < _expected_n:
                            # Verify re-expansion gives correct count
                            try:
                                _p2 = _CP.from_str(
                                    _rescued_cif, occupancy_tolerance=100.0)
                                _s2 = _p2.parse_structures(primitive=False)
                                if _s2 and abs(len(_s2[0]) - _expected_n) <= 1:
                                    print(f"  CifWriter rescue succeeded: "
                                          f"SG {_rsg}, {len(_rsites)} sites "
                                          f"(symprec={_symprec})", flush=True)
                                    # Keep CifWriter's cell and positions
                                    # together; patching only the cell axes
                                    # recreates the mixed-setting CIF bug.
                                    return _rescued_cif
                            except BaseException:
                                pass
                    except BaseException:
                        continue
                print(f"  CifWriter rescue failed — will use synthetic CIF "
                      f"(may have wrong atom positions)", flush=True)
        except BaseException as _e:
            print(f"  from_spacegroup validation error: {_e}", flush=True)

    # ── Safety check: detect if "reduction" actually returned full-cell sites ──
    # If the reduction fell through all fallbacks and returned raw CIF sites
    # that are actually the full unit cell (no symops in the original CIF,
    # e.g. from a pymatgen Structure dump), writing those sites with the
    # declared SG > 1 would cause GSAS-II to double-expand.
    #
    # Detection: compare site count with what CifParser(primitive=False) gives.
    # If they match AND the CIF lacks symmetry operations, the sites are full-cell.
    # In this case, use the ORIGINAL CIF directly — it already has the correct
    # atom positions and (if from CifWriter) proper symmetry operations.
    if sites and sg > 1 and cif_text:
        has_symops = ('_symmetry_equiv_pos_as_xyz' in cif_text or
                      '_space_group_symop_operation_xyz' in cif_text or
                      '_space_group_symop.operation_xyz' in cif_text)
        if not has_symops:
            # The original CIF has no symmetry operations — check if sites
            # are the full cell by expanding the original CIF
            try:
                from pymatgen.io.cif import CifParser
                try:
                    _parser = CifParser.from_str(cif_text, occupancy_tolerance=100.0)
                except (AttributeError, TypeError):
                    import tempfile as _tf, os as _os
                    _fd, _tmp = _tf.mkstemp(suffix='.cif')
                    with _os.fdopen(_fd, 'w') as _f:
                        _f.write(cif_text)
                    _parser = CifParser(_tmp, occupancy_tolerance=100.0)
                    _os.unlink(_tmp)
                _structs = _parser.parse_structures(primitive=False)
                if _structs:
                    full_n = len(_structs[0])
                    if len(sites) >= full_n:
                        # Sites ARE the full cell — DO NOT write them with SG > 1.
                        # Instead, try to use CifWriter to generate a proper CIF
                        # with asymmetric unit + symmetry operations.
                        print(f"  _build_conventional_cif: detected full-cell "
                              f"sites ({len(sites)}) with no symops for "
                              f"SG {sg}. Attempting CifWriter rescue...",
                              flush=True)
                        try:
                            from pymatgen.io.cif import CifWriter
                            # Use the already-parsed structure
                            struct0 = _structs[0]
                            for symprec in (0.01, 0.05, 0.1, 0.2):
                                try:
                                    writer = CifWriter(struct0, symprec=symprec)
                                    import tempfile as _tf2, os as _os2
                                    _fd2, _tmp2 = _tf2.mkstemp(suffix='.cif')
                                    with _os2.fdopen(_fd2, 'w') as _f2:
                                        pass
                                    writer.write_file(_tmp2)
                                    with open(_tmp2, 'r') as _f2:
                                        rescued_cif = _f2.read()
                                    _os2.unlink(_tmp2)
                                    # Validate: parse the rescued CIF
                                    rescued_parsed = parse_cif(rescued_cif)
                                    rescued_sg = rescued_parsed.get(
                                        'spacegroup_number', 1)
                                    rescued_sites = rescued_parsed.get(
                                        'sites') or []
                                    if (rescued_sites
                                            and len(rescued_sites) < full_n
                                            and rescued_sg > 1):
                                        print(f"  CifWriter rescue succeeded: "
                                              f"SG {rescued_sg}, "
                                              f"{len(rescued_sites)} asym sites "
                                              f"(symprec={symprec})", flush=True)
                                        # Keep CifWriter's cell and positions
                                        # together. If CifWriter standardized
                                        # axes, patching only the cell lengths
                                        # makes an inconsistent refinement CIF.
                                        return rescued_cif
                                except BaseException:
                                    continue
                        except Exception as e_cw:
                            print(f"  CifWriter rescue failed: {e_cw}",
                                  flush=True)

                        # CifWriter rescue failed — last resort: write with P1
                        # so GSAS-II doesn't expand. This will give correct
                        # peak positions but GSAS-II won't apply symmetry
                        # constraints on cell parameters or atoms.
                        print(f"  LAST RESORT: writing full-cell sites with "
                              f"P1 to prevent double-expansion. GSAS-II will "
                              f"not apply symmetry constraints.", flush=True)
                        sg = 1
                        hm = 'P 1'
            except BaseException:
                pass

    # ── Diagnostic: print what GSAS-II will actually receive ──────────
    print(f"  _build_conventional_cif [{formula or '?'}]:", flush=True)
    print(f"    Cell: a={a:.4f} b={b:.4f} c={c:.4f} "
          f"α={al:.2f} β={be:.2f} γ={ga:.2f}", flush=True)
    print(f"    SG: {sg} ({hm})", flush=True)
    print(f"    Asymmetric sites: {len(sites) if sites else 0}", flush=True)
    if sites:
        for _el, _x, _y, _z, _occ in sites:
            print(f"      {_el:3s}  ({_x:.6f}, {_y:.6f}, {_z:.6f})  "
                  f"occ={_occ:.4f}", flush=True)

    lines = [
        'data_phase',
        f"_chemical_formula_sum '{formula}'" if formula else '',
        f"_cell_formula_units_Z {Z}" if Z else '',
        f'_cell_length_a {a:.5f}',
        f'_cell_length_b {b:.5f}',
        f'_cell_length_c {c:.5f}',
        f'_cell_angle_alpha {al:.3f}',
        f'_cell_angle_beta {be:.3f}',
        f'_cell_angle_gamma {ga:.3f}',
        f'_symmetry_Int_Tables_number {sg}',
        f"_symmetry_space_group_name_H-M '{hm}'",
        f"_space_group_IT_number {sg}",
        f"_space_group_name_H-M_alt '{hm}'",
    ]

    if sites:
        lines += [
            '',
            'loop_',
            '_atom_site_label',
            '_atom_site_type_symbol',
            '_atom_site_fract_x',
            '_atom_site_fract_y',
            '_atom_site_fract_z',
            '_atom_site_occupancy',
        ]
        site_counts = {}
        for el, x, y, z, occ in sites:
            # Generate unique labels like W1, W2, O1, O2...
            site_counts[el] = site_counts.get(el, 0) + 1
            label = f'{el}{site_counts[el]}'
            # Clamp negative zeros to 0.0 — GSAS-II CIF reader rejects "-0.000000"
            x = abs(x) if abs(x) < 1e-8 else x
            y = abs(y) if abs(y) < 1e-8 else y
            z = abs(z) if abs(z) < 1e-8 else z
            lines.append(f'{label}  {el}  {x:.6f}  {y:.6f}  {z:.6f}  {occ:.4f}')

    return '\n'.join(ln for ln in lines if ln is not None) + '\n'


# ─────────────────────────────────────────────────────────────────────────────
# NON-NEGATIVE LEAST SQUARES (NNLS)
# ─────────────────────────────────────────────────────────────────────────────

def _nnls(A, b):
    """Solve  min‖Ax − b‖²  subject to  x ≥ 0.

    Uses scipy.optimize.nnls if available, otherwise falls back to the
    Lawson-Hanson active-set algorithm (pure numpy).
    """
    try:
        from scipy.optimize import nnls as _scipy_nnls
        x, _ = _scipy_nnls(A, b)
        return x
    except ImportError:
        pass
    # ── Fallback: Lawson-Hanson active-set NNLS ──────────────────────
    m, n = A.shape
    P = set()               # passive set (variables free to be > 0)
    Z = set(range(n))       # zero set (variables fixed at 0)
    x = np.zeros(n)
    max_outer = 3 * n
    for _ in range(max_outer):
        w = A.T @ (b - A @ x)  # negative gradient
        if not Z or max(w[j] for j in Z) <= 1e-10:
            break
        # Move variable with largest positive gradient into passive set
        t = max(Z, key=lambda j: w[j])
        P.add(t); Z.discard(t)
        # Inner loop: solve unconstrained on P, handle negatives
        for __ in range(max_outer):
            cols = sorted(P)
            s = np.zeros(n)
            s[cols] = np.linalg.lstsq(A[:, cols], b, rcond=None)[0]
            if all(s[j] >= 0 for j in cols):
                x = s
                break
            # Find the tightest step that hits zero
            alpha = 1.0
            j_remove = None
            for j in cols:
                if s[j] < 0 and x[j] > x[j] - s[j]:
                    a = x[j] / (x[j] - s[j])
                    if a < alpha:
                        alpha = a
                        j_remove = j
            if j_remove is None:
                x = s
                break
            x = x + alpha * (s - x)
            x[j_remove] = 0.0
            Z.add(j_remove); P.discard(j_remove)
        else:
            x = np.maximum(x, 0.0)  # safety clamp
    return x


# ─────────────────────────────────────────────────────────────────────────────
# PER-PHASE PATTERN FROM GSAS-II REFLECTIONS
# ─────────────────────────────────────────────────────────────────────────────

def _compute_raw_phase_profile(tt_arr, refs, U_deg, V_deg, W_deg,
                                X_deg, Y_deg, gaussian_only=False):
    """Compute a raw (unscaled) profile shape for one phase.

    Returns an array of the same length as *tt_arr* whose values are
    proportional to the diffracted intensity at each 2θ point.  The
    absolute scale is arbitrary — what matters is the *relative* shape
    compared to other phases so that the proportional decomposition
    can correctly partition the total pattern.

    Parameters
    ----------
    gaussian_only : bool
        If True, use pure Gaussian profiles (eta=0) with a tighter
        window (3×FWHM).  This eliminates the long Lorentzian tails
        that cause cross-phase leakage artifacts in the proportional
        decomposition, while preserving accurate peak localisation
        for genuinely overlapping reflections.

    refs : list of [two_theta, d, (h,k,l), intensity_weight] from
           generate_reflections (with structure-factor filtering).
    """
    n = len(tt_arr)
    pattern = np.zeros(n, dtype=np.float64)

    use_tch = (X_deg != 0.0 or Y_deg != 0.0)

    for tt_pos, d, hkl, weight in refs:
        if weight <= 0:
            continue

        if use_tch:
            fwhm, eta = tch_fwhm_eta(tt_pos, U_deg, V_deg, W_deg,
                                      X_deg, Y_deg)
        else:
            fwhm = max(caglioti_fwhm(tt_pos, U_deg, V_deg, W_deg), 0.005)
            eta = 0.5

        if gaussian_only:
            eta = 0.0

        sigma_g = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        gamma_l = fwhm / 2.0
        # Gaussian-only: 3×FWHM is sufficient (G ≈ 10⁻¹¹ at edge).
        # Pseudo-Voigt: 6×FWHM captures Lorentzian tails.
        window = 3.0 * fwhm if gaussian_only else 6.0 * fwhm

        msk = np.abs(tt_arr - tt_pos) < window
        if not msk.any():
            continue

        dx = tt_arr[msk] - tt_pos
        G = np.exp(-0.5 * (dx / sigma_g) ** 2)
        if eta > 0:
            L = 1.0 / (1.0 + (dx / gamma_l) ** 2)
            prof = eta * L + (1.0 - eta) * G
        else:
            prof = G
        # Normalise
        area = np.trapz(prof, tt_arr[msk]) if np.trapz(prof, tt_arr[msk]) > 0 else 1.0
        pattern[msk] += weight * prof / area

    return pattern


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _auto_select_bg_coeffs(tt, y_obs, phases, wavelength, tt_min, tt_max,
                           user_n=6):
    """Automatically select the number of Chebyshev background coefficients.

    Philosophy: use enough terms to accurately model the instrument
    background, but not so many that the polynomial develops oscillations
    (humps/dips) that absorb peak tails or noise.  Even clean crystalline
    samples need ~6 terms to capture instrument effects (air scatter,
    Bremsstrahlung, sample holder scattering, beam spillover at low 2θ).
    Going below 6 starves the background of flexibility and degrades Rwp.

    Strategy:
      1. Identify "peak-free" regions by masking out zones within ±1.0° of
         expected reflection positions (wider exclusion to avoid peak tails
         biasing the background estimate).
      2. Fit peak-free intensities with Chebyshev polynomials of order 6, 10.
      3. Default to 6 terms for clean crystalline patterns.
      4. Step up to 10 only if there is a clear amorphous hump (DW < 0.5
         AND broad feature detected AND >20% RMS improvement from 6→10),
         as seen in supported catalysts (WC/SiO₂, Pt/Al₂O₃).
      5. Only go higher (14, 20) if the user explicitly requests it.

    Returns the recommended n_bg_coeffs (int).
    """
    from .crystallography import generate_reflections

    # Build a mask of peak-free zones
    mask_range = (tt >= tt_min) & (tt <= tt_max)
    tt_r = tt[mask_range]
    y_r = y_obs[mask_range]

    if len(tt_r) < 20 or tt_max <= tt_min:
        return user_n

    # Find all expected peak positions from ALL phases
    peak_positions = []
    for ph in phases:
        try:
            sys_ = (ph.get('system') or 'triclinic').lower()
            sg = ph.get('spacegroup_number', 1)
            refs = generate_reflections(
                ph.get('a', 4), ph.get('b', ph.get('a', 4)),
                ph.get('c', ph.get('a', 4)),
                ph.get('alpha', 90), ph.get('beta', 90), ph.get('gamma', 90),
                sys_, sg, wavelength, tt_min, tt_max, hkl_max=10)
            peak_positions.extend([r[0] for r in refs])
        except Exception:
            pass

    # Mask: True = peak-free (at least 1.0° away from any reflection)
    # Using wider exclusion (1.0° vs 0.5°) to keep peak tails out of the
    # baseline estimate — Lorentzian tails extend far from peak centres
    # and can bias the polynomial if included.
    peak_free = np.ones(len(tt_r), dtype=bool)
    for pos in peak_positions:
        peak_free &= np.abs(tt_r - pos) > 1.0

    # Need enough peak-free points to fit
    n_free = int(np.sum(peak_free))
    if n_free < 15:
        print(f"  Auto-BG: only {n_free} peak-free points — using "
              f"n_bg={user_n}", flush=True)
        return user_n

    tt_free = tt_r[peak_free]
    y_free = y_r[peak_free]

    # Normalize 2theta to [-1, 1] for Chebyshev stability
    x_free = 2.0 * (tt_free - tt_min) / (tt_max - tt_min) - 1.0

    # Fit with 3 terms (constant + linear + quadratic — handles most
    # instrument-induced slopes like beam spillover at low angles)
    try:
        coeffs_3 = np.polynomial.chebyshev.chebfit(x_free, y_free, 2)
        fitted_3 = np.polynomial.chebyshev.chebval(x_free, coeffs_3)
        resid_3 = y_free - fitted_3
        rms_3 = float(np.sqrt(np.mean(resid_3 ** 2)))
    except Exception:
        print(f"  Auto-BG: Chebyshev fit failed — using n_bg=6", flush=True)
        return 6

    # Durbin-Watson statistic: tests for autocorrelation in residuals.
    # DW ≈ 2.0 → no autocorrelation (residuals are random noise → good fit).
    # DW < 1.0 → strong positive autocorrelation (systematic misfit →
    #             the background has curvature that 3 terms can't capture).
    dw = 2.0
    if len(resid_3) > 2:
        diff_sq = float(np.sum(np.diff(resid_3) ** 2))
        sum_sq = float(np.sum(resid_3 ** 2))
        if sum_sq > 0:
            dw = diff_sq / sum_sq

    # Also check if the baseline has a broad hump: fit with 6 terms and
    # see if any higher-order coefficients are large relative to the
    # constant term (indicating genuine curvature, not noise).
    has_broad_feature = False
    rms_6 = rms_3
    try:
        coeffs_6 = np.polynomial.chebyshev.chebfit(x_free, y_free, 5)
        fitted_6 = np.polynomial.chebyshev.chebval(x_free, coeffs_6)
        resid_6 = y_free - fitted_6
        rms_6 = float(np.sqrt(np.mean(resid_6 ** 2)))
        # Check if higher-order coefficients are significant
        if abs(coeffs_6[0]) > 0:
            max_higher = max(abs(c) for c in coeffs_6[3:])
            if max_higher > abs(coeffs_6[0]) * 0.05:
                has_broad_feature = True
    except Exception:
        pass

    # Also fit with 10 terms to check if amorphous hump needs even more
    # flexibility (e.g. SiO₂, Al₂O₃ support backgrounds)
    rms_10 = rms_6
    try:
        coeffs_10 = np.polynomial.chebyshev.chebfit(x_free, y_free, 9)
        fitted_10 = np.polynomial.chebyshev.chebval(x_free, coeffs_10)
        resid_10 = y_free - fitted_10
        rms_10 = float(np.sqrt(np.mean(resid_10 ** 2)))
    except Exception:
        pass

    # Decision logic — 6 terms is the floor for all samples.
    # Even clean crystalline patterns need 6 Chebyshev terms to capture
    # instrument background effects (air scatter, Bremsstrahlung, sample
    # holder scattering, beam spillover at low 2θ).  Dropping to 3 starves
    # the background and degrades Rwp.
    #
    # Tier 1 (6 terms): clean crystalline baseline — standard default
    #   → No broad features, or mild curvature captured by 6 terms
    #
    # Tier 2 (10 terms): significant background structure
    #   Triggered by ANY of these conditions:
    #   (a) DW < 0.1 on the 3-term residuals — extreme autocorrelation
    #       indicates the background has curvature that 6 terms likely
    #       cannot capture either.  A DW this low means the background
    #       model is catastrophically underfitting.
    #   (b) DW < 0.5 AND broad feature AND >15% RMS improvement from
    #       6→10 terms (relaxed from 20% — borderline cases like
    #       WC + W2C patterns with DW=0.04 were being missed).
    #   (c) DW(6) < 0.3 on the 6-term residuals — even after doubling
    #       the terms, there's still strong autocorrelation.
    #
    # Beyond 10: user must explicitly select (14, 20) from the dropdown
    improvement_6_to_10 = (rms_6 - rms_10) / max(rms_6, 1e-6)

    # Also compute DW on the 6-term residuals to check if 6 is adequate
    dw_6 = 2.0
    try:
        if len(resid_6) > 2:
            _diff_sq_6 = float(np.sum(np.diff(resid_6) ** 2))
            _sum_sq_6 = float(np.sum(resid_6 ** 2))
            if _sum_sq_6 > 0:
                dw_6 = _diff_sq_6 / _sum_sq_6
    except Exception:
        pass

    # Default to 6 Chebyshev terms.  Escalate to 10 when the DW/RMS
    # tests indicate the background is underfit.
    best_n = 6
    _escalation_reason = None

    # Tier 2 triggers — escalate from 6 to 10 only when BOTH
    # autocorrelation AND RMS improvement confirm the background is
    # genuinely underfit.  DW alone is not sufficient because
    # nanocrystalline patterns with closely spaced broad peaks have
    # small peak-free zones where residuals are artificially
    # autocorrelated from peak tails, not from background underfitting.
    if dw < 0.1 and has_broad_feature and improvement_6_to_10 > 0.10:
        # (a) Extreme autocorrelation + broad feature + 10% RMS improvement
        best_n = 10
        _escalation_reason = (f"DW(3)={dw:.3f} < 0.1 + broad_feature + "
                              f"improv={improvement_6_to_10:.1%} > 10%")
    elif dw < 0.5 and has_broad_feature and improvement_6_to_10 > 0.20:
        # (b) Moderate autocorrelation + broad feature + strong improvement
        best_n = 10
        _escalation_reason = (f"DW(3)={dw:.3f} < 0.5 + broad_feature + "
                              f"improv={improvement_6_to_10:.1%} > 20%")
    elif dw_6 < 0.3 and improvement_6_to_10 > 0.15:
        # (c) Strong autocorrelation on 6-term residuals + improvement
        best_n = 10
        _escalation_reason = (f"DW(6)={dw_6:.3f} < 0.3 + "
                              f"improv={improvement_6_to_10:.1%} > 15%")

    # User override: if user explicitly requested more than auto-selected
    if user_n > best_n:
        best_n = user_n
        _escalation_reason = f"user requested {user_n}"

    print(f"  Auto-BG: peak-free analysis — RMS(6)={rms_6:.2f}, "
          f"RMS(10)={rms_10:.2f}, DW(3)={dw:.2f}, DW(6)={dw_6:.2f}, "
          f"broad_feature={has_broad_feature}, "
          f"improv_6→10={improvement_6_to_10:.1%} → n_bg={best_n}"
          + (f" [{_escalation_reason}]" if _escalation_reason else ""),
          flush=True)
    return best_n


def _write_xye(path, tt, y_obs, sigma):
    """Write a .xye file (2theta  intensity  sigma) for GSAS-II import."""
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        for i in range(len(tt)):
            f.write(f"{tt[i]:.6f}  {y_obs[i]:.4f}  {sigma[i]:.4f}\n")


def _estimate_profile_params(tt, y_obs):
    """Estimate initial Caglioti U, V, W from observed peak widths.

    Finds the strongest peaks in the data, measures their approximate FWHM,
    and fits the Caglioti equation FWHM² = U tan²θ + V tanθ + W to get
    reasonable starting values.  Returns (U, V, W) in centidegrees² / centideg²
    (GSAS-II internal units), or the module defaults if estimation fails.

    This is critical for complex phases like W2C where the default U/V/W
    may be far from the true broadening, causing GSAS-II to converge to
    a local minimum.
    """
    try:
        from scipy.signal import find_peaks

        # Subtract a simple baseline (percentile)
        baseline = np.percentile(y_obs, 10)
        y_corr = y_obs - baseline

        # Find prominent peaks (height > 10% of max, well separated)
        height_thresh = 0.1 * y_corr.max()
        step = float(tt[1] - tt[0]) if len(tt) > 1 else 0.02
        distance = max(1, int(0.5 / step))  # at least 0.5° apart
        peaks, props = find_peaks(y_corr, height=height_thresh, distance=distance)

        if len(peaks) < 2:
            return DEFAULT_U, DEFAULT_V, DEFAULT_W

        # Measure FWHM for each peak
        fwhm_data = []  # (tan_theta, fwhm_deg)
        for pk in peaks:
            half_max = y_corr[pk] / 2.0
            # Walk left
            left = pk
            while left > 0 and y_corr[left] > half_max:
                left -= 1
            # Walk right
            right = pk
            while right < len(y_corr) - 1 and y_corr[right] > half_max:
                right += 1
            fwhm_deg = float(tt[right] - tt[left])
            if 0.02 < fwhm_deg < 3.0:  # reasonable range
                theta_rad = math.radians(float(tt[pk]) / 2.0)
                tan_th = math.tan(theta_rad)
                fwhm_data.append((tan_th, fwhm_deg))

        if len(fwhm_data) < 2:
            return DEFAULT_U, DEFAULT_V, DEFAULT_W

        # Fit Caglioti: FWHM² = U tan²θ + V tanθ + W
        # In GSAS-II, U/V/W are in centideg² so FWHM is in centideg.
        # We work in degrees then convert.
        tan_arr = np.array([d[0] for d in fwhm_data])
        fwhm_sq = np.array([d[1] ** 2 for d in fwhm_data])

        # Build design matrix [tan²θ, tanθ, 1]
        A = np.column_stack([tan_arr ** 2, tan_arr, np.ones_like(tan_arr)])
        try:
            result, _, _, _ = np.linalg.lstsq(A, fwhm_sq, rcond=None)
            U_deg2, V_deg2, W_deg2 = float(result[0]), float(result[1]), float(result[2])
        except np.linalg.LinAlgError:
            return DEFAULT_U, DEFAULT_V, DEFAULT_W

        # Convert degrees² → centideg²  (multiply by 10000)
        U_cdeg2 = U_deg2 * 10000.0
        V_cdeg2 = V_deg2 * 10000.0
        W_cdeg2 = W_deg2 * 10000.0

        # Sanity clamp: U and W should be positive, V typically negative
        U_cdeg2 = max(0.1, min(U_cdeg2, 500.0))
        V_cdeg2 = max(-200.0, min(V_cdeg2, 200.0))
        W_cdeg2 = max(0.1, min(W_cdeg2, 500.0))

        return U_cdeg2, V_cdeg2, W_cdeg2

    except Exception:
        return DEFAULT_U, DEFAULT_V, DEFAULT_W


def _estimate_lorentzian_params(tt, y_obs, U_cdeg2, V_cdeg2, W_cdeg2):
    """Estimate initial Lorentzian profile parameters X and Y from peak shapes.

    Measures the Lorentzian fraction of observed peaks by comparing their
    half-max width to their quarter-max width (Gaussian ratio ~1.48,
    Lorentzian ratio ~1.73).  Decomposes each peak's FWHM into Gaussian and
    Lorentzian components using the Thompson-Cox-Hastings (TCH) relationship,
    then fits H_L = X*tan(theta) + Y/cos(theta) across peaks.

    Returns (X, Y) in centidegrees, or moderate defaults if estimation fails.
    Critical for carbides/oxides with significant crystallite-size (Y) and
    micro-strain (X) Lorentzian broadening.
    """
    try:
        from scipy.signal import find_peaks

        baseline = np.percentile(y_obs, 10)
        y_corr = y_obs - baseline

        height_thresh = 0.1 * y_corr.max()
        step = float(tt[1] - tt[0]) if len(tt) > 1 else 0.02
        distance = max(1, int(0.5 / step))
        peaks, _ = find_peaks(y_corr, height=height_thresh, distance=distance)

        if len(peaks) < 2:
            return DEFAULT_X, DEFAULT_Y

        lor_data = []  # (tan_theta, cos_theta, H_L_deg)
        for pk in peaks:
            half_max = y_corr[pk] / 2.0
            quarter_max = y_corr[pk] / 4.0

            # Measure half-max width
            left_h, right_h = pk, pk
            while left_h > 0 and y_corr[left_h] > half_max:
                left_h -= 1
            while right_h < len(y_corr) - 1 and y_corr[right_h] > half_max:
                right_h += 1
            fwhm_deg = float(tt[right_h] - tt[left_h])

            # Measure quarter-max width
            left_q, right_q = pk, pk
            while left_q > 0 and y_corr[left_q] > quarter_max:
                left_q -= 1
            while right_q < len(y_corr) - 1 and y_corr[right_q] > quarter_max:
                right_q += 1
            fwqm_deg = float(tt[right_q] - tt[left_q])

            if fwhm_deg < 0.02 or fwhm_deg > 3.0 or fwqm_deg <= fwhm_deg:
                continue

            # Estimate Lorentzian fraction eta from FWQM/FWHM ratio.
            # Gaussian: FWQM/FWHM = sqrt(2) ≈ 1.414
            # Lorentzian: FWQM/FWHM = sqrt(3) ≈ 1.732
            # Linear interpolation gives eta (Lorentzian fraction).
            ratio = fwqm_deg / fwhm_deg
            # Linear interpolation for eta
            eta = (ratio - 1.4142) / (1.7321 - 1.4142)
            eta = max(0.0, min(eta, 1.0))

            if eta < 0.01:
                continue  # essentially pure Gaussian, no Lorentzian info

            # TCH decomposition: total FWHM ≈ eta*H_L + (1-eta)*H_G (simplified)
            # More precise: use the 5th-order TCH formula inverse.
            # For a practical estimate, the pseudo-Voigt approximation gives:
            #   H_L ≈ FWHM * eta (leading-order)
            #   H_G ≈ FWHM * (1 - eta)
            # Subtract the Gaussian Caglioti contribution to isolate Lorentzian.
            theta_rad = math.radians(float(tt[pk]) / 2.0)
            tan_th = math.tan(theta_rad)
            cos_th = math.cos(theta_rad)

            # Caglioti Gaussian FWHM² in deg² (convert from centideg²)
            H_G_sq_deg2 = (U_cdeg2 * tan_th**2 + V_cdeg2 * tan_th
                           + W_cdeg2) / 10000.0
            H_G_deg = math.sqrt(max(H_G_sq_deg2, 0.001))

            # Lorentzian FWHM: from total FWHM and Gaussian contribution
            # Using TCH relation: H_total^5 ≈ H_G^5 + ... + H_L^5
            # Simplified: H_L ≈ max(0, FWHM - H_G) when eta > 0
            # More robust: H_L = FWHM * eta (pseudo-Voigt definition)
            H_L_deg = fwhm_deg * eta
            if H_L_deg > 0.005:
                lor_data.append((tan_th, cos_th, H_L_deg))

        if len(lor_data) < 2:
            return DEFAULT_X, DEFAULT_Y

        # Fit H_L = X*tan(theta) + Y/cos(theta)
        # X and Y are in degrees here; convert to centideg at the end.
        tan_arr = np.array([d[0] for d in lor_data])
        inv_cos_arr = np.array([1.0 / d[1] for d in lor_data])
        hl_arr = np.array([d[2] for d in lor_data])

        A = np.column_stack([tan_arr, inv_cos_arr])
        try:
            result, _, _, _ = np.linalg.lstsq(A, hl_arr, rcond=None)
            X_deg, Y_deg = float(result[0]), float(result[1])
        except np.linalg.LinAlgError:
            return DEFAULT_X, DEFAULT_Y

        # Convert degrees → centidegrees (multiply by 100)
        X_cdeg = X_deg * 100.0
        Y_cdeg = Y_deg * 100.0

        # Sanity clamp: both should be non-negative, moderate magnitude
        X_cdeg = max(0.0, min(X_cdeg, 50.0))
        Y_cdeg = max(0.0, min(Y_cdeg, 50.0))

        return X_cdeg, Y_cdeg

    except Exception:
        return DEFAULT_X, DEFAULT_Y


# Cu Kα1/Kα2 doublet constants (NIST values).  Used by _write_instprm
# whenever the incoming wavelength is close to the Cu Kα weighted mean.
# Most Cu-source diffractometers export either the weighted average
# (~1.5418 Å) or Kα1-only (~1.5406 Å); both fall inside the window below.
# If the user truly wants a monochromated single-wavelength fit they can
# pass `kalpha2=False` explicitly to _write_instprm.
CU_KALPHA1_A     = 1.540593   # Å  (Cu Kα1, NIST)
CU_KALPHA2_A     = 1.544414   # Å  (Cu Kα2, NIST)
CU_KALPHA2_RATIO = 0.5        # I(Kα2) / I(Kα1), canonical for Cu
CU_DETECT_MIN_A  = 1.537      # auto-detect window: any wavelength
CU_DETECT_MAX_A  = 1.548      # in [1.537, 1.548] Å is treated as Cu.


def _is_cu_kalpha(wavelength):
    """True if *wavelength* (Å) is within the Cu Kα auto-detect window."""
    try:
        w = float(wavelength)
    except (TypeError, ValueError):
        return False
    return CU_DETECT_MIN_A <= w <= CU_DETECT_MAX_A


def _write_instprm(work_dir, wavelength, polariz=None, sh_l=None,
                   u=None, v=None, w=None, x=None, y=None,
                   kalpha2=None, zero_seed=None):
    """Write a minimal GSAS-II .instprm file. Returns path.

    Uses module-level DEFAULT_* constants unless overridden.
    U, V, W, X, Y are initial guesses that GSAS-II will refine.
    Polariz. and SH/L are NOT refined — they should match the instrument.
    zero_seed: starting value for Zero parameter (default 0.0).

    Wavelength handling
    -------------------
    - When *kalpha2* is True, write a Cu Kα1/Kα2 doublet profile
      (`Lam1`, `Lam2`, `I(L2)/I(L1)`).  The passed-in *wavelength*
      argument is ignored in favour of the NIST Cu Kα constants; this
      matches what a Kβ-filtered Cu lab diffractometer actually emits.
    - When *kalpha2* is False, write a single-wavelength profile (`Lam`).
    - When *kalpha2* is None (default), auto-detect: if *wavelength*
      falls within the Cu Kα window (~1.54 Å) write the doublet;
      otherwise write single wavelength.  This fixes the silent misfit
      for unmonochromated Cu data without requiring a UI change.
    """
    path = os.path.join(work_dir, 'instrument.instprm')
    pol = polariz if polariz is not None else DEFAULT_POLARIZ
    shl = sh_l if sh_l is not None else DEFAULT_SH_L
    _u = u if u is not None else DEFAULT_U
    _v = v if v is not None else DEFAULT_V
    _w = w if w is not None else DEFAULT_W
    _x = x if x is not None else DEFAULT_X
    _y = y if y is not None else DEFAULT_Y

    # Decide wavelength mode.
    if kalpha2 is None:
        use_doublet = _is_cu_kalpha(wavelength)
    else:
        use_doublet = bool(kalpha2)

    if use_doublet:
        # Cu Kα1 / Kα2 doublet.  GSAS-II treats Type:PXC with Lam1/Lam2
        # keys as a two-wavelength calculation and generates every peak
        # as a Kα1/Kα2 pair at the correct Bragg-law-scaled 2θ positions.
        wavelength_lines = [
            f'Lam1:{CU_KALPHA1_A:.6f}',
            f'Lam2:{CU_KALPHA2_A:.6f}',
            f'I(L2)/I(L1):{CU_KALPHA2_RATIO:.4f}',
        ]
        print(f"  instprm: Cu Kα1/Kα2 doublet "
              f"(Lam1={CU_KALPHA1_A:.6f}, Lam2={CU_KALPHA2_A:.6f}, "
              f"I2/I1={CU_KALPHA2_RATIO}).", flush=True)
    else:
        wavelength_lines = [f'Lam:{wavelength:.6f}']
        print(f"  instprm: single wavelength "
              f"(Lam={wavelength:.6f}).", flush=True)

    # Zero seed: starting value for GSAS-II's Zero parameter refinement.
    # This is just an initial guess — Zero is still refined freely in
    # Stage 2 and Stage 6.  The value is now instrument-profile-specific:
    #   synergy_s: -0.25° (from Synergy-Dualflex .par)
    #   smartlab:   0.0° (no prior calibration data)
    _ZERO_SEED = zero_seed if zero_seed is not None else 0.0
    lines = [
        '#GSAS-II instrument parameter file; do not add/delete items!',
        'Type:PXC',
    ] + wavelength_lines + [
        f'Zero:{_ZERO_SEED}',
        f'Polariz.:{pol}',
        f'U:{_u}',
        f'V:{_v}',
        f'W:{_w}',
        f'X:{_x}',
        f'Y:{_y}',
        'Z:0.0',
        f'SH/L:{shl}',
        'Azimuth:0.0',
    ]
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(lines) + '\n')
    return path


def _write_temp_cif(cif_text, phase_name='phase', work_dir=None, index=0):
    """Write CIF text to a temporary file. Returns path.

    If *work_dir* is given the file is created there (avoids Windows
    short-path / permission issues with the system temp directory).
    The *index* parameter ensures unique filenames when multiple phases
    share the same name.
    """
    # GSAS-II's CIF validator opens files through the Windows default text
    # codec in some installs (cp1252), so UTF-8-only punctuation in comments
    # can make an otherwise valid CIF unreadable. Keep the file handed to
    # GSAS-II strictly ASCII; crystallographic tags/values are ASCII anyway.
    _ascii_map = {
        '\u2192': '->', '\u2013': '-', '\u2014': '-',
        '\u2018': "'", '\u2019': "'", '\u201c': '"', '\u201d': '"',
        '\u03b1': 'alpha', '\u03b2': 'beta', '\u03b3': 'gamma',
        '\u00b5': 'u', '\u00b0': ' deg', '\u00b2': '^2',
        '\u00b3': '^3', '\u00c5': 'A',
    }
    cif_text = ''.join(_ascii_map.get(ch, ch) for ch in str(cif_text))
    cif_text = cif_text.encode('ascii', 'replace').decode('ascii')

    if work_dir is not None:
        # Sanitise phase_name for use as a filename component
        safe = "".join(c if c.isalnum() or c in ('_', '-') else '_'
                       for c in phase_name)
        path = os.path.join(work_dir, f'gsas_{index}_{safe}.cif')
        with open(path, 'w', encoding='ascii', newline='\n') as f:
            f.write(cif_text)
        return path
    # Fallback: system temp
    fd, path = tempfile.mkstemp(suffix='.cif', prefix=f'gsas_{index}_{phase_name}_')
    with os.fdopen(fd, 'w', encoding='ascii', newline='\n') as f:
        f.write(cif_text)
    return path


def _extract_profile_params(phase_obj):
    """
    Extract TCH profile parameters from a GSAS-II phase object.
    GSAS-II stores profile coefficients in the histogram data.
    Returns dict with crystallite size and microstrain.

    GSAS-II Size data structure (isotropic):
      Size = ['isotropic', [size_value, refine_flag], ...]
    where size_value is in MICRONS (µm), NOT Angstroms.
    Convert: size_A = size_um * 10000

    Mustrain data structure (isotropic):
      Mustrain = ['isotropic', [strain_value, refine_flag], ...]
    where strain_value is in units of 10^-6 (micro-strain).
    """
    try:
        # GSAS-II stores profile params per histogram in the phase
        hapData = list(phase_obj.data['Histograms'].values())[0]
        size_data = hapData.get('Size', [])
        strain_data = hapData.get('Mustrain', [])

        cryst_size_A = None
        if size_data and len(size_data) > 1:
            size_val = float(size_data[1][0]) if size_data[1][0] > 0 else None
            if size_val is not None and size_val > 0:
                # GSAS-II stores Size in µm; convert to Å (×10000)
                cryst_size_A = size_val * 10000.0
                # Sanity check: crystallite size < 5 Å is physically
                # impossible (smaller than a unit cell).  This catches
                # cases where Size was never refined and GSAS-II returned
                # its unphysical initial value.
                if cryst_size_A < 5.0:
                    print(f"  WARNING: GSAS-II crystallite size = "
                          f"{cryst_size_A:.1f} Å ({size_val:.6f} µm) — "
                          f"unphysically small, discarding (will use "
                          f"Y-parameter fallback).", flush=True)
                    cryst_size_A = None

        microstrain = None
        if strain_data and len(strain_data) > 1:
            microstrain = float(strain_data[1][0])

        return {
            'crystallite_size_A': cryst_size_A,
            'microstrain_pct': microstrain,
        }
    except Exception:
        return {'crystallite_size_A': None, 'microstrain_pct': None}


def _extract_instrument_params(histogram):
    """Extract instrument/profile parameters from a GSAS-II histogram."""
    try:
        inst_params = histogram.data['Instrument Parameters'][0]
        U = float(inst_params.get('U', [0, 0])[1])
        V = float(inst_params.get('V', [0, 0])[1])
        W = float(inst_params.get('W', [0, 0])[1])
        X = float(inst_params.get('X', [0, 0])[1])
        Y = float(inst_params.get('Y', [0, 0])[1])
        return {'U': U, 'V': V, 'W': W, 'X': X, 'Y': Y}
    except Exception:
        return {'U': 0, 'V': 0, 'W': 0.1, 'X': 0, 'Y': 0}


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def run_gsas2(tt, y_obs, sigma, phases, wavelength,
              tt_min=None, tt_max=None, n_bg_coeffs=6,
              max_cycles=10, progress_callback=None,
              instprm_file=None, polariz=None, sh_l=None,
              auto_bg=True, seed_params=None, options=None,
              instrument=None, instrument_reason=None):
    """
    Run GSAS-II Rietveld refinement via GSASIIscriptable.

    Parameters
    ----------
    tt, y_obs, sigma : np.ndarray — data arrays
    phases   : list of dicts — must include 'cif_text' for atom positions
    wavelength : float (Å)
    tt_min/max : float — fitting range
    n_bg_coeffs : int — number of background coefficients (starting point
        when auto_bg=True, exact value when auto_bg=False)
    max_cycles : int — max refinement cycles
    instprm_file : str, optional — path to user-provided .instprm file;
        if given, this is used INSTEAD of all provisional instrument
        defaults.  This is the single highest-priority setting — a
        measured standard .instprm pins U/V/W/X/Y/SH/L to physical
        values and makes sample parameters uniquely determinable.
    polariz : float, optional — monochromator polarization (default 0.99)
    sh_l : float, optional — Finger-Cox-Jephcoat asymmetry (default 0.002)
    auto_bg : bool — if True, automatically determine the optimal number
        of background coefficients from the data (default True)
    options : dict, optional — backend refinement options with safe
        defaults.  Callers can ignore this entirely; the backend will
        behave predictably with the defaults below.

        Supported keys:
          geometry : str — "capillary" (default) or "bragg_brentano".
              Controls which displacement correction is refined:
              capillary → DisplaceY, bragg_brentano → DisplaceX.
          preferred_orientation : str — "auto" (default), "off", or
              "force".  "auto" enables March-Dollase for hexagonal/
              trigonal phases only.  "off" disables for all.  "force"
              enables for all phases.
          refine_xyz : bool — whether to refine atom positions (Stage 5).
              Default False.  Without a measured .instprm, XYZ
              refinement can absorb errors from polarization,
              absorption, PO, background, or wrong phase choice.
          background_mode : str — "auto" (default), "fixed6", "fixed10",
              or "fixed14".  "auto" uses DW/RMS analysis to select
              between 6 and 10.  Fixed modes bypass auto-selection.
          exclude_regions : list of (lo, hi) tuples — 2θ ranges to
              mask out before fitting (detector artifacts, holder
              peaks, low-angle junk).
          phase_sensitivity : bool — if True, run diagnostic variants
              (PO on/off, XYZ on/off, BG 6/10) and report stability.
              Default False.  Slow but useful for automated screening.

    Returns dict compatible with Le Bail / Rietveld output (same keys).
    """
    print("\n=== GSAS-II REFINEMENT START ===", flush=True)
    if not _GSASII_AVAILABLE:
        raise RuntimeError(
            f"GSAS-II is not installed or not importable in this Python environment. "
            f"Install with: git clone https://github.com/AdvancedPhotonSource/GSAS-II.git "
            f"&& cd GSAS-II && pip install .\n"
            f"Import error: {_GSASII_IMPORT_ERROR}")

    # ── Resolve instrument profile ─────────────────────────────────────────
    # The instrument profile bundles geometry, displacement model, polariz,
    # SH/L, PO default, zero seed, and sigma inflation into one coherent
    # set of defaults.  Explicit params (polariz=, sh_l=) override profile.
    instrument = instrument or DEFAULT_INSTRUMENT
    instrument_reason = instrument_reason or 'default'
    profile = INSTRUMENT_PROFILES.get(instrument, INSTRUMENT_PROFILES[DEFAULT_INSTRUMENT])

    print(f"  Instrument: {profile['label']}", flush=True)
    print(f"    Selection: {instrument_reason}", flush=True)
    print(f"    Displacement model: {profile['displacement_param']}",
          flush=True)
    # Note: "Profile:" line is printed AFTER auto-locate below,
    # so it correctly reflects whether a measured .instprm was found.

    # ── Parse options with safe defaults ───────────────────────────────────
    # Options can override profile defaults where needed.
    options = options or {}
    geometry             = options.get('geometry', profile['geometry'])
    preferred_orientation = options.get('preferred_orientation',
                                        profile['preferred_orientation_default'])
    refine_xyz           = options.get('refine_xyz', False)
    background_mode      = options.get('background_mode', 'auto')
    exclude_regions      = options.get('exclude_regions', [])
    phase_sensitivity    = options.get('phase_sensitivity', False)
    # Verification mode: skip Stage 3 (cell), Stage 4 (Uiso+MD), Stage 4b
    # (size), and disable cell refinement in Stage 6.  First-pass test for
    # a new sample/CIF combination — refines only background, scales,
    # displacement, and Y, so position/width problems become visible
    # without being absorbed by the more flexible parameters.
    verification_mode    = options.get('verification_mode', False)
    # When verification_mode AND verify_refine_cell are both True, Stage 3
    # (per-phase cell with divergence checks) still skips, but Stage 6
    # gets cell refinement back on top of bg + scales + displacement + Y.
    # Use this once verification_mode has shown peaks land in the right
    # place and you want to let the cell shrink/relax.  Uiso, size, MD,
    # and XYZ remain off.
    verify_refine_cell   = options.get('verify_refine_cell', False)
    # Phase isolation: per-phase patterns from GSAS-II ycalc instead of
    # manual reconstruction.  Default ON — internally consistent with
    # GSAS-II RefList ticks (see tick-source switch in per-phase loop).
    phase_isolation_opt  = bool(options.get('phase_isolation', True))
    # Refine March-Dollase preferred orientation for hexagonal/trigonal
    # phases during verification mode.  Without this flag, PO is set up
    # at ratio=1.0 (no texture) but never refines, so hex platelet
    # texture (e.g. WC along [001]) cannot be captured.  When True,
    # forces preferred_orientation to 'auto' (hex/trig only) and lets
    # Stage 6 refine the MD ratio alongside cell.
    verify_refine_po_opt = bool(options.get('verify_refine_po', False))
    if verify_refine_po_opt and verification_mode:
        if preferred_orientation == 'off':
            preferred_orientation = 'auto'
            print(f"  PO override: preferred_orientation forced to 'auto' "
                  f"because verify_refine_po=True.", flush=True)
    # Swap position-handle: free Zero, fix DisplaceX/Y at zero.  Use when
    # the offset is a real diffractometer/wavelength miscalibration that
    # the measured-instprm Zero cannot capture for this sample (e.g. when
    # DisplaceY refuses to move from 0 because the position residual
    # belongs to Zero, not to sample displacement).  Overrides the
    # measured-instprm rule that locks Zero.
    use_zero_not_displace_opt = bool(
        options.get('verify_use_zero_not_displace', False))
    if use_zero_not_displace_opt:
        print(f"  Position handle: refining Zero, fixing DisplaceX/Y = 0 "
              f"(verify_use_zero_not_displace=True).", flush=True)
    # Constrain W2C cell to uniform-volume scaling (Branch B).
    # Stage 6 refines (a, b, c) freely; afterward, post-process W2C
    # cell to enforce (a, b, c) = (s·a₀, s·b₀, s·c₀) where
    # s = (V_refined/V_start)^(1/3).  Tests whether the unconstrained
    # anisotropic refinement was capturing real strain or just exploiting
    # the data-resolution-limited a/b correlation.  Compare χ² against
    # the unconstrained run: similar χ² → uniform scaling is enough,
    # the anisotropy was noise; much worse χ² → real anisotropic strain.
    uniform_cell_w2c_opt = bool(
        options.get('verify_cell_uniform_w2c', False))
    if uniform_cell_w2c_opt:
        print(f"  Cell constraint: W2C will be post-scaled to uniform "
              f"contraction after Stage 6 "
              f"(verify_cell_uniform_w2c=True).", flush=True)
    # Diagnostic: free Lorentzian strain X in Stage 2 and Stage 6.
    # Normally X is pinned to the measured-instprm value (Si standard
    # has X = 0).  When this is True, X is added to the refine list
    # alongside Y (and Zero, if zero-swap is also on).  Use to test
    # whether residual Lorentzian peak shape needs strain-like broadening.
    refine_x_opt = bool(options.get('verify_refine_x', False))
    if refine_x_opt:
        print(f"  Diagnostic: X (Lorentzian strain) will refine in "
              f"Stages 2 and 6 (verify_refine_x=True).", flush=True)
    # Y controls — orthogonal to X and Zero handles.  Fix-Y dominates
    # refine-Y when both are set.  Non-negative constraint applies
    # post-Stage-6 only when Y refines.
    fix_y_opt           = bool(options.get('verify_fix_y', False))
    y_fixed_value_opt   = options.get('verify_y_fixed_value', None)
    y_nonneg_opt        = bool(options.get('verify_y_nonnegative', False))
    if fix_y_opt:
        print(f"  Y control: fixing Y at requested value "
              f"(verify_fix_y=True, value={y_fixed_value_opt}).",
              flush=True)
    if y_nonneg_opt:
        print(f"  Y control: post-Stage-6 will clamp Y to >= 0 "
              f"(verify_y_nonnegative=True).", flush=True)
    _y_nonnegative_clamped = False
    # Explicit toggles for structural refinement (Phase B).  These
    # OVERRIDE the verification_mode default of "off" — they let the
    # user run a near-verification recipe with one extra structural
    # parameter active, without leaving verification_mode entirely.
    refine_uiso_opt     = bool(options.get('verify_refine_uiso', False))
    refine_size_opt     = bool(options.get('verify_refine_size', False))
    if refine_uiso_opt:
        print(f"  Structural: Uiso (atomic displacement) will refine in "
              f"Stage 4 (verify_refine_uiso=True).", flush=True)
    if refine_size_opt:
        print(f"  Structural: HAP Size will refine in Stage 4b "
              f"(verify_refine_size=True).", flush=True)
    # Tick source: default Python phase_refs (filtered by intensity),
    # opt-in to GSAS-II RefList (all reflections including weak ones).
    use_gsas_ticks_opt  = bool(options.get('use_gsas_ref_ticks', False))
    # Fix WC PO: keep March-Dollase enabled at a held value instead of
    # refining.  Use when adding Uiso to a recipe where MD was refining
    # — Uiso × MD correlation is the main SVD source in those runs.
    fix_po_opt          = bool(options.get('verify_fix_po', False))
    try:
        po_fixed_value_opt = float(
            options.get('verify_po_fixed_value', 0.905) or 0.905)
    except (TypeError, ValueError):
        po_fixed_value_opt = 0.905
    if fix_po_opt:
        print(f"  PO control: WC PO held at MD ratio = "
              f"{po_fixed_value_opt:.4f} (verify_fix_po=True).",
              flush=True)
    # Per-phase refinement options (dynamic, phase-agnostic).
    # phase_options is a list (one entry per phase, indexed) of dicts:
    #   {"refine_size": bool, "refine_mustrain": bool,
    #    "po_mode": "off"|"fixed"|"refined",
    #    "po_value": float, "po_axis": [h, k, l] or None}
    # When provided, the backend applies these per-phase by index instead
    # of the legacy name-based WC/W2C options.
    phase_options_raw = options.get('phase_options', None)
    if phase_options_raw is None:
        phase_options_list = []
    elif isinstance(phase_options_raw, list):
        phase_options_list = phase_options_raw
    elif isinstance(phase_options_raw, dict):
        # Convert {"0": {...}, "1": {...}} → list
        try:
            _max_idx = max(int(k) for k in phase_options_raw.keys())
            phase_options_list = [
                phase_options_raw.get(str(i), {})
                for i in range(_max_idx + 1)]
        except Exception:
            phase_options_list = []
    else:
        phase_options_list = []
    if phase_options_list:
        print(f"  Per-phase options received for {len(phase_options_list)} "
              f"phase(s):", flush=True)
        for _pi, _popts in enumerate(phase_options_list):
            print(f"    Phase {_pi}: {_popts}", flush=True)

    # Legacy phase-specific HAP options — kept as a fallback for callers
    # that haven't migrated to phase_options.  Detection by case-insensitive
    # name match: any phase whose name contains "wc" but not "w2c" → WC;
    # contains "w2c" → W2C.
    refine_wc_size_opt      = bool(options.get('verify_refine_wc_size', False))
    refine_w2c_size_opt     = bool(options.get('verify_refine_w2c_size', False))
    refine_wc_mustrain_opt  = bool(options.get('verify_refine_wc_mustrain', False))
    refine_w2c_mustrain_opt = bool(options.get('verify_refine_w2c_mustrain', False))
    if any([refine_wc_size_opt, refine_w2c_size_opt,
            refine_wc_mustrain_opt, refine_w2c_mustrain_opt]):
        print(f"  HAP per-phase (legacy name-based): "
              f"WC size={refine_wc_size_opt}, "
              f"W2C size={refine_w2c_size_opt}, "
              f"WC mustrain={refine_wc_mustrain_opt}, "
              f"W2C mustrain={refine_w2c_mustrain_opt}",
              flush=True)

    def _phase_opts_for(idx, name):
        """Return effective per-phase options dict.  phase_options_list
        wins over legacy WC/W2C name-based options.  Returns a dict with
        all keys present (None/defaults for unspecified)."""
        defaults = {
            'refine_size': False,
            'refine_mustrain': False,
            'po_mode': 'off',
            'po_value': 1.0,
            'po_axis': None,
            'uniform_cell': False,
        }
        # Layer 1: legacy name-based fallback
        _n = (name or '').lower()
        if 'w2c' in _n:
            if refine_w2c_size_opt:
                defaults['refine_size'] = True
            if refine_w2c_mustrain_opt:
                defaults['refine_mustrain'] = True
        elif 'wc' in _n:
            if refine_wc_size_opt:
                defaults['refine_size'] = True
            if refine_wc_mustrain_opt:
                defaults['refine_mustrain'] = True
        # Layer 2: explicit phase_options entry overrides
        if idx < len(phase_options_list):
            _explicit = phase_options_list[idx] or {}
            if isinstance(_explicit, dict):
                for k in defaults:
                    if k in _explicit:
                        defaults[k] = _explicit[k]
        return defaults

    _uniform_cell_phase_idxs = set()
    for _idx_uc, _ph_uc in enumerate(phases):
        _name_uc = (_ph_uc.get('name') or _ph_uc.get('formula') or '')
        if bool(_phase_opts_for(_idx_uc, _name_uc).get('uniform_cell')):
            _uniform_cell_phase_idxs.add(_idx_uc)
    if _uniform_cell_phase_idxs:
        uniform_cell_w2c_opt = True
        print(f"  Cell constraint: uniform-cell scaling requested for "
              f"phase index(es) {sorted(_uniform_cell_phase_idxs)}.",
              flush=True)
    if verification_mode:
        if verify_refine_cell:
            print("  VERIFICATION + CELL MODE: Uiso, MD, size, and XYZ "
                  "still skipped.  Cell refines in Stage 6 only.  "
                  "Refining bg + scales + displacement + Y + cell.",
                  flush=True)
        else:
            print("  VERIFICATION MODE: cell, Uiso, MD, size, and XYZ "
                  "stages will be skipped.  Refining bg + scales + "
                  "displacement + Y only.  Add cell refinement in a "
                  "follow-up run once peaks land in the right place.",
                  flush=True)
        refine_xyz = False  # never refine atom positions in verification mode

    # Profile-derived defaults (explicit params take priority)
    if polariz is None:
        polariz = profile['polariz']
    if sh_l is None:
        sh_l = profile['sh_l']

    # Sigma inflation factor — profile-specific
    _SIGMA_INFLATION_FACTOR = profile.get('sigma_inflation_K', 1.0)

    print(f"  Options: geometry={geometry}, PO={preferred_orientation}, "
          f"XYZ={refine_xyz}, BG={background_mode}, "
          f"exclude={len(exclude_regions)} regions, "
          f"sigma_K={_SIGMA_INFLATION_FACTOR}", flush=True)

    # Validate: all phases need CIF text
    missing = [ph.get('name', '?') for ph in phases if not ph.get('cif_text')]
    if missing:
        raise ValueError(
            f"GSAS-II refinement requires CIF with atom coordinates for all phases. "
            f"Missing for: {', '.join(missing)}.")

    if tt_min is None: tt_min = float(tt.min())
    if tt_max is None: tt_max = float(tt.max())

    # Apply range mask
    mask = (tt >= tt_min) & (tt <= tt_max)
    tt_r = tt[mask]
    y_r = y_obs[mask]
    sig_r = sigma[mask] if sigma is not None else np.sqrt(np.maximum(y_r, 1.0))

    # ── Apply exclusion regions ─────────────────────────────────────────
    # Mask out user-specified 2θ ranges (detector artifacts, holder peaks,
    # low-angle junk, or known non-phase features) before writing XYE.
    if exclude_regions:
        keep = np.ones(len(tt_r), dtype=bool)
        for lo, hi in exclude_regions:
            keep &= ~((tt_r >= lo) & (tt_r <= hi))
        n_excluded = int(np.sum(~keep))
        tt_r = tt_r[keep]
        y_r = y_r[keep]
        sig_r = sig_r[keep]
        print(f"  Exclusion regions: masked {n_excluded} points in "
              f"{len(exclude_regions)} region(s).", flush=True)

    if progress_callback:
        progress_callback('GSAS-II: setting up project...')

    # ── Create temporary files ───────────────────────────────────────────
    # Use a temp dir inside the app directory to avoid spaces / short-path
    # issues in the user-profile temp folder (GSAS-II can't read files
    # whose path contains spaces).
    _app_tmp = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            os.pardir, os.pardir, '.gsas_tmp')
    _app_tmp = os.path.normpath(_app_tmp)
    os.makedirs(_app_tmp, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix='gsas2_', dir=_app_tmp)
    gpx_path = os.path.join(work_dir, 'refine.gpx')
    data_path = os.path.join(work_dir, 'data.xye')
    # Use user-provided .instprm if available; otherwise generate defaults.
    #
    # Seed-handling note (Phase I v5): the in-house Rietveld seed pipeline
    # fits peak widths using a SINGLE wavelength model, so its U/V/W/X/Y
    # already have the unresolved Kα1/Kα2 doublet broadening baked in.
    # When GSAS-II is about to model the doublet EXPLICITLY, those seeds
    # are ~2× too wide — using them causes GSAS-II to refine U downward
    # and cascades into scale/cell errors (observed on the first Cu
    # doublet run: wt% flipped from 19/81 to 2/98).  Solution: ignore
    # seed_params whenever the doublet will be written, and let the
    # observation-based estimator (_estimate_profile_params) derive seed
    # widths from the actual data.  The estimator already measures the
    # total apparent FWHM, which GSAS-II will then deconvolve correctly
    # once the Kα1/Kα2 pair is being forward-modelled.
    _use_doublet = _is_cu_kalpha(wavelength)
    _seed_active = (seed_params is not None) and (not _use_doublet)
    if _seed_active is False and seed_params is not None and _use_doublet:
        print(f"  Skipping single-wavelength seed params "
              f"(Cu Kα1/Kα2 doublet active — seeds would be "
              f"~2× too wide for the doublet model).", flush=True)

    # Auto-locate measured .instprm from instrument profile if not
    # explicitly provided.  Searches relative to the toolkit root
    # (two levels up from this module file).
    if not instprm_file:
        _profile_instprm = profile.get('instprm_filename')
        if _profile_instprm:
            _module_dir = os.path.dirname(os.path.abspath(__file__))
            _toolkit_root = os.path.dirname(os.path.dirname(_module_dir))
            _candidate = os.path.join(_toolkit_root, _profile_instprm)
            if os.path.isfile(_candidate):
                instprm_file = _candidate
                print(f"  Auto-located measured .instprm: {_candidate}",
                      flush=True)

    # Report profile status AFTER auto-locate, so it's accurate
    if instprm_file and os.path.isfile(instprm_file):
        print(f"    Profile: measured .instprm (overrides provisional "
              f"defaults)", flush=True)
    else:
        print(f"    Profile: temporary fallback, not measured .instprm",
              flush=True)

    _measured_instprm = False  # True when U/V/W/X/Y should be FIXED

    if instprm_file and os.path.isfile(instprm_file):
        import shutil
        instprm_path = os.path.join(work_dir, 'instrument.instprm')
        shutil.copy2(instprm_file, instprm_path)
        _measured_instprm = True
        print(f"Using MEASURED instrument parameters: {instprm_file}",
              flush=True)
        print(f"  → U/V/W/X/SH/L and Zero fixed; Y may refine for "
              f"sample size broadening."
              f"\n  → Refining: displacement, Y, scale, cell, BG"
              f"  (and Size/Uiso unless verification_mode).",
              flush=True)
        # Populate est_* with defaults for downstream reset targets; the
        # user .instprm controls actual refinement starting values.
        est_u, est_v, est_w = DEFAULT_U, DEFAULT_V, DEFAULT_W
        est_x, est_y = DEFAULT_X, DEFAULT_Y
    elif _seed_active:
        # Use pre-refinement seed values (e.g. from in-house Rietveld).
        # Already in GSAS-II units (centideg² for U/V/W, centideg for X/Y).
        est_u = float(seed_params.get('U', DEFAULT_U))
        est_v = float(seed_params.get('V', DEFAULT_V))
        est_w = float(seed_params.get('W', DEFAULT_W))
        est_x = float(seed_params.get('X', DEFAULT_X))
        est_y = float(seed_params.get('Y', DEFAULT_Y))
        # Sanity clamp — wider than _estimate_profile_params because the
        # in-house Rietveld can legitimately find large U/W for
        # nanocrystalline materials (e.g. U=3.5 deg² = 35000 centideg²
        # for ~9 nm crystallites).
        est_u = max(0.1, min(est_u, 50000.0))      # centideg²
        est_v = max(-10000.0, min(est_v, 10000.0))
        est_w = max(0.1, min(est_w, 50000.0))
        est_x = max(0.0, min(est_x, 500.0))        # centideg
        est_y = max(0.0, min(est_y, 500.0))
        print(f"  Using seeded profile params: U={est_u:.2f}, "
              f"V={est_v:.2f}, W={est_w:.2f}, X={est_x:.2f}, "
              f"Y={est_y:.2f}", flush=True)
        instprm_path = _write_instprm(work_dir, wavelength,
                                       polariz=polariz, sh_l=sh_l,
                                       u=est_u, v=est_v, w=est_w,
                                       x=est_x, y=est_y,
                                       zero_seed=profile.get('zero_seed', 0.0))
    else:
        # Estimate initial profile parameters from observed peak widths.
        # This path is used when no user instprm and no usable seeds —
        # either because seed_params is None, or because the Cu doublet
        # is active (see _seed_active logic above).
        est_u, est_v, est_w = _estimate_profile_params(tt_r, y_r)
        est_x, est_y = _estimate_lorentzian_params(tt_r, y_r,
                                                    est_u, est_v, est_w)
        instprm_path = _write_instprm(work_dir, wavelength,
                                       polariz=polariz, sh_l=sh_l,
                                       u=est_u, v=est_v, w=est_w,
                                       x=est_x, y=est_y,
                                       zero_seed=profile.get('zero_seed', 0.0))
    # Save initial profile estimates for use as reset targets if
    # parameters diverge during refinement.  Better than resetting
    # to arbitrary small constants (e.g. 0.01) which are far from
    # any physical value and corrupt subsequent stages.
    _init_W = est_w if not (instprm_file and os.path.isfile(instprm_file)) else DEFAULT_W
    _init_U = est_u if not (instprm_file and os.path.isfile(instprm_file)) else DEFAULT_U
    _write_xye(data_path, tt_r, y_r, sig_r)

    cif_paths = []
    gsas_cif_texts = []  # Store CIF text sent to GSAS-II for each phase
    for i, ph in enumerate(phases):
        # Build a synthetic CIF from the phase dict's (conventional) cell
        # parameters.  This guarantees GSAS-II sees the correct space group
        # and cell geometry even when the original CIF used a primitive
        # setting (common with Materials Project data).
        cif_for_gsas = _build_conventional_cif(ph)
        gsas_cif_texts.append(cif_for_gsas)

        # ── W2C / mp-2034 emergency guard ─────────────────────────────
        # The fixture pipeline must produce a Pbcn (SG 60) CIF for
        # mp-2034.  If anything in app.py / cif_cache / _build_conventional_cif
        # has rewritten the CIF into P1/full-cell, refining further is
        # silently meaningless.  Stop the run with a clear error
        # instead of producing garbage results.
        _phase_id = str(ph.get('cod_id', ph.get('mp_id', '')))
        if _phase_id == 'mp-2034':
            _has_pbcn = ('Pbcn' in cif_for_gsas)
            _has_sg60 = ('_symmetry_Int_Tables_number   60' in cif_for_gsas
                         or '_symmetry_Int_Tables_number 60' in cif_for_gsas
                         or '_space_group_IT_number 60' in cif_for_gsas
                         or '_space_group_IT_number   60' in cif_for_gsas)
            if not (_has_pbcn and _has_sg60):
                raise RuntimeError(
                    f"W2C mp-2034 is not using the Pbcn fixture CIF.  "
                    f"Expected fixtures/w2c_pbcn_mp_2034.cif before "
                    f"GSAS-II refinement, but the CIF being written "
                    f"to {ph.get('name', 'phase')!r} contains "
                    f"Pbcn={_has_pbcn}, SG60={_has_sg60}.  Refinement "
                    f"aborted to prevent silent garbage output.  "
                    f"Check app.py phase_hint merge logic and "
                    f"cif_cache.cached_fetch_mp fixture priority."
                )
            print(f"  ✓ W2C mp-2034 fixture guard passed: Pbcn + SG 60.",
                  flush=True)

        cif_path = _write_temp_cif(cif_for_gsas, ph.get('name', 'phase'),
                                   work_dir=work_dir, index=i)
        cif_paths.append(cif_path)
        # ── Diagnostic: dump CIF content for debugging ────────────────
        print(f"\n  === CIF for phase {i} ({ph.get('name', '?')}) ===",
              flush=True)
        print(f"  File: {cif_path}", flush=True)
        # Print first 30 lines of CIF (enough to see cell + atoms)
        _cif_lines = cif_for_gsas.strip().split('\n')
        for _ln in _cif_lines[:40]:
            print(f"    {_ln}", flush=True)
        if len(_cif_lines) > 40:
            print(f"    ... ({len(_cif_lines) - 40} more lines)", flush=True)
        print(f"  === end CIF ===\n", flush=True)

    try:
        # ── Build GSAS-II project ────────────────────────────────────────
        try:
            gpx = G2sc.G2Project(newgpx=gpx_path)
        except Exception as e:
            raise RuntimeError(
                f"GSAS-II failed to create a project. This usually means GSAS-II's "
                f"compiled Fortran extensions are not installed (requires a Fortran "
                f"compiler during installation).\n\n"
                f"Try using the built-in Rietveld refinement instead, or install "
                f"GSAS-II via: conda install gsas2pkg -c briantoby\n\n"
                f"Original error: {e}"
            ) from e

        # Add histogram (powder data)
        histogram = gpx.add_powder_histogram(
            data_path, instprm_path,
            databank=None,
        )

        # Set data range
        histogram.data['Limits'] = [[tt_min, tt_max], [tt_min, tt_max]]

        # ── Sigma inflation for 2D-integrated / capillary data ─────────
        # Data exported from CrysAlisPro powder generation (or any
        # similar 2D→1D ring-integration pipeline: Dioptas, pyFAI, GSAS-II
        # 2D, etc.) carries per-point σ that reflects ONLY counting
        # statistics across thousands of detector pixels per 2θ bin.
        # Those σ are typically 20–40× tighter than classical √I Poisson
        # on the summed ring — statistically correct for random noise,
        # but they do not include experimental systematics (beam drift,
        # capillary transparency tail, sample-height micromotion, detector
        # response nonuniformities).
        #
        # GSAS-II weights each point by w = 1/σ², so tight σ inflate
        # reduced χ² and GoF even when the structural model is sound.
        # For the WC/W2C Synergy-Dualflex dataset, reduced χ² ≈ 87 while
        # Rwp ≈ 4.8% — a classic signature of underestimated σ, not a
        # broken fit.  Inflate σ by a fixed factor via wtFactor (which
        # multiplies w uniformly: setting wtFactor = 1/K² inflates σ by
        # K and divides reduced χ² by K²).  Rwp is invariant under
        # uniform σ scaling because it is a ratio.
        #
        # Sigma inflation factor K is now instrument-profile-specific:
        #   synergy_s:  K=5.0 (2D-integrated σ underestimates uncertainty)
        #   smartlab:   K=1.0 (σ ≈ √I is already correct for BB)
        # The value was resolved from the instrument profile at the top
        # of this function into _SIGMA_INFLATION_FACTOR.
        try:
            _wt = 1.0 / (_SIGMA_INFLATION_FACTOR ** 2)
            histogram.data['wtFactor'] = _wt
            print(f"  σ inflation: wtFactor = {_wt:.4f} "
                  f"(σ × {_SIGMA_INFLATION_FACTOR:.1f}, reduced χ² ÷ "
                  f"{_SIGMA_INFLATION_FACTOR ** 2:.0f}). "
                  f"Rwp is invariant under this rescaling.", flush=True)
        except (KeyError, TypeError) as _wt_e:
            print(f"  WARNING: could not set wtFactor — {_wt_e}",
                  flush=True)

        # Fix histogram scale to 1.0 — NEVER refine it.
        # GSAS-II has N+1 scale parameters (N phase scales + 1 histogram
        # scale).  Only N are independent.  Refining all N+1 creates a
        # perfect correlation (100%) that causes SVD singularity and the
        # refinement gets stuck with zero peak intensity.
        # Standard practice: fix histogram scale, refine only phase scales.
        try:
            histogram.data['Sample Parameters']['Scale'] = [1.0, False]
        except (KeyError, TypeError):
            pass

        # Set background
        bkg_data = histogram.data['Background']
        bg_init = float(np.percentile(y_r, 2))
        bkg_data[0] = ['chebyschev-1', True, n_bg_coeffs,
                        bg_init] + [0.0] * (n_bg_coeffs - 1)

        # Add phases from CIF files — always embed the space group number in
        # the GSAS-II phasename so that two phases with the same formula
        # (e.g. W Pm-3n vs W Im-3m) are never treated as the same phase.
        gsas_phases = []
        used_names = set()
        for i, (ph, cif_path) in enumerate(zip(phases, cif_paths)):
            formula  = ph.get('formula', '') or ph.get('name', '') or 'Phase'
            sg_num   = ph.get('spacegroup_number', 1)
            # Always include SG number for unambiguous internal naming
            phasename = f"{formula}_sg{sg_num}"
            # Deduplicate in case two phases have the same formula + SG (rare)
            if phasename in used_names:
                phasename = f"{phasename}_{i+1}"
            used_names.add(phasename)
            try:
                phase_obj = gpx.add_phase(
                    cif_path,
                    phasename=phasename,
                    histograms=[histogram],
                )
            except Exception as e1:
                # Synthetic CIF failed — fall back to original CIF text
                orig_cif = ph.get('cif_text', '')
                if orig_cif:
                    try:
                        orig_path = _write_temp_cif(
                            orig_cif, ph.get('name', 'phase'),
                            work_dir=work_dir, index=100+i)
                        cif_paths.append(orig_path)  # for cleanup
                        phase_obj = gpx.add_phase(
                            orig_path,
                            phasename=phasename,
                            histograms=[histogram],
                        )
                        warnings.warn(
                            f"Synthetic CIF failed for '{ph.get('name', '?')}', "
                            f"using original CIF text (error: {e1})")
                    except Exception as e2:
                        cif_size = os.path.getsize(cif_path) if os.path.exists(cif_path) else -1
                        cif_head = ''
                        if os.path.exists(cif_path):
                            with open(cif_path, 'r', encoding='utf-8') as _f:
                                cif_head = _f.read(500)
                        raise RuntimeError(
                            f"GSAS-II could not read CIF for phase '{ph.get('name', '?')}' "
                            f"(file: {cif_path}, size: {cif_size} bytes).\n"
                            f"CIF preview:\n{cif_head}\n\n"
                            f"Original error: {e2}"
                        ) from e2
                else:
                    cif_size = os.path.getsize(cif_path) if os.path.exists(cif_path) else -1
                    cif_head = ''
                    if os.path.exists(cif_path):
                        with open(cif_path, 'r', encoding='utf-8') as _f:
                            cif_head = _f.read(500)
                    raise RuntimeError(
                        f"GSAS-II could not read CIF for phase '{ph.get('name', '?')}' "
                        f"(file: {cif_path}, size: {cif_size} bytes).\n"
                        f"CIF preview:\n{cif_head}\n\n"
                        f"Original error: {e1}"
                    ) from e1
            gsas_phases.append(phase_obj)

        # ── Diagnostic: print GSAS-II's expanded atoms for each phase ─────
        for _pi, _po in enumerate(gsas_phases):
            _atoms = _po.data.get('Atoms', [])
            _gcell = _po.data['General']['Cell']
            print(f"  GSAS-II phase {_pi} ({_po.name}): "
                  f"{len(_atoms)} atoms, SG={_po.data['General'].get('SGData', {}).get('SpGrp', '?')}, "
                  f"cell=[{_gcell[1]:.4f}, {_gcell[2]:.4f}, {_gcell[3]:.4f}, "
                  f"{_gcell[4]:.2f}, {_gcell[5]:.2f}, {_gcell[6]:.2f}]",
                  flush=True)
            for _ai, _at in enumerate(_atoms[:20]):  # print first 20 atoms
                _lbl = _at[0] if len(_at) > 0 else '?'
                _el  = _at[1] if len(_at) > 1 else '?'
                _x   = _at[3] if len(_at) > 3 else '?'
                _y   = _at[4] if len(_at) > 4 else '?'
                _z   = _at[5] if len(_at) > 5 else '?'
                print(f"    atom {_ai}: {_lbl} ({_el}) at "
                      f"({_x}, {_y}, {_z})", flush=True)

        # ── CIF import validation gate ─────────────────────────────────────
        # Compare what GSAS-II actually imported against what the search-row
        # phase metadata SAID it should be.  Catches the W2C-as-P1 class of
        # bug (where CIF construction silently produces a structure GSAS-II
        # interprets as a different space group / setting / atom count).
        # This is non-fatal — we log a clear warning and continue, leaving
        # the user to abort if the mismatch matters for their fit.
        _validation_warnings = []
        _validation_failures = []
        for _pi, (_po, _ph_input) in enumerate(zip(gsas_phases, phases)):
            _expected_sg = int(_ph_input.get('spacegroup_number', 0) or 0)
            _expected_a  = float(_ph_input.get('a', 0) or 0)
            _expected_b  = float(_ph_input.get('b', _expected_a) or _expected_a)
            _expected_c  = float(_ph_input.get('c', _expected_a) or _expected_a)
            _name        = _ph_input.get('formula') or _ph_input.get('name') or f'phase{_pi}'

            try:
                _imported_sgdata = _po.data['General'].get('SGData', {}) or {}
                _imported_sg_str = (_imported_sgdata.get('SpGrp')
                                    or _imported_sgdata.get('SpcGrp')
                                    or '?').strip()
                _imported_sg_num = _hm_symbol_to_number(_imported_sg_str)
                _imported_atoms = len(_po.data.get('Atoms', []))
                _imported_cell  = _po.data['General']['Cell']
                _ia = float(_imported_cell[1])
                _ib = float(_imported_cell[2])
                _ic = float(_imported_cell[3])
                _expected_atoms = None
                try:
                    _gsas_cif_text = (gsas_cif_texts[_pi]
                                      if _pi < len(gsas_cif_texts) else '')
                    _expected_atoms = len(
                        parse_cif(_gsas_cif_text).get('sites') or [])
                except Exception:
                    _expected_atoms = None
            except Exception as _e:
                print(f"  ⚠ Validation [{_pi} {_name}]: could not read "
                      f"GSAS-II state: {_e}", flush=True)
                continue

            # Space group sanity check.  P1 (1) when something else was
            # expected → mostly-fatal: refinement is now fitting a
            # different model than the user thinks.
            if _expected_sg > 1 and _imported_sg_num == 1:
                _msg = (f"phase {_pi} ({_name}): expected SG "
                        f"{_expected_sg} but GSAS-II imported as "
                        f"P1 (SG 1).  CIF round-trip likely lost "
                        f"symmetry — refinement model is NOT what "
                        f"the search row indicated.")
                _validation_warnings.append(_msg)
                _validation_failures.append(_msg)
                print(f"  ⚠ VALIDATION FAIL: {_msg}", flush=True)
            elif _expected_sg > 1 and _imported_sg_num != _expected_sg:
                _msg = (f"phase {_pi} ({_name}): expected SG "
                        f"{_expected_sg}, GSAS-II imported SG "
                        f"{_imported_sg_num} ('{_imported_sg_str}').  "
                        f"Possible setting transformation (CifWriter "
                        f"may have standardized axes).")
                _validation_warnings.append(_msg)
                print(f"  ⚠ VALIDATION WARN: {_msg}", flush=True)

            # Cell sanity: each axis within ±5% of expected, or a
            # PERMUTATION of (a, b, c) within ±5%.  Permutation is
            # legal for orthorhombic CifWriter standardization but
            # worth flagging.
            if (_expected_atoms is not None and _expected_atoms > 0
                    and _imported_atoms != _expected_atoms):
                _msg = (f"phase {_pi} ({_name}): expected "
                        f"{_expected_atoms} CIF atom site(s), but GSAS-II "
                        f"imported {_imported_atoms}.")
                if _expected_sg > 1 and _imported_atoms > _expected_atoms:
                    _validation_failures.append(_msg)
                    print(f"  VALIDATION FAIL: {_msg}", flush=True)
                else:
                    _validation_warnings.append(_msg)
                    print(f"  VALIDATION WARN: {_msg}", flush=True)

            def _close(x, y, tol=0.05):
                if x <= 0 or y <= 0:
                    return False
                return abs(x - y) / max(x, y) <= tol
            if _expected_a > 0 and _expected_b > 0 and _expected_c > 0:
                _direct = (_close(_ia, _expected_a)
                           and _close(_ib, _expected_b)
                           and _close(_ic, _expected_c))
                _expected_sorted = sorted([_expected_a, _expected_b, _expected_c])
                _imported_sorted = sorted([_ia, _ib, _ic])
                _permuted = all(_close(x, y) for x, y in
                                zip(_expected_sorted, _imported_sorted))
                if not _direct and _permuted:
                    _msg = (f"phase {_pi} ({_name}): cell axes match "
                            f"under permutation but not direct.  "
                            f"Expected ({_expected_a:.4f}, "
                            f"{_expected_b:.4f}, {_expected_c:.4f}) → "
                            f"got ({_ia:.4f}, {_ib:.4f}, {_ic:.4f}).  "
                            f"OK for orthorhombic if positions were "
                            f"transformed together, but verify the "
                            f"refined cell maps to your intended axes.")
                    _validation_warnings.append(_msg)
                    print(f"  ⚠ VALIDATION WARN: {_msg}", flush=True)
                elif not _direct and not _permuted:
                    _msg = (f"phase {_pi} ({_name}): cell mismatch.  "
                            f"Expected ({_expected_a:.4f}, "
                            f"{_expected_b:.4f}, {_expected_c:.4f}), "
                            f"got ({_ia:.4f}, {_ib:.4f}, {_ic:.4f}).")
                    _validation_warnings.append(_msg)
                    _validation_failures.append(_msg)
                    print(f"  ⚠ VALIDATION FAIL: {_msg}", flush=True)

        if _validation_failures:
            detail = "\n".join(f"  - {m}" for m in _validation_failures)
            raise RuntimeError(
                "GSAS-II CIF import validation failed before refinement.\n"
                f"{detail}\n"
                "The phase model imported by GSAS-II does not match the "
                "selected CIF/metadata, so refinement was aborted.")

        if _validation_warnings:
            print(f"  ⚠ {len(_validation_warnings)} validation issue(s) "
                  f"detected before refinement.  Continuing anyway — "
                  f"check the warnings above and the refined cell vs. "
                  f"expected at the end.", flush=True)
        else:
            print(f"  ✓ CIF import validation: all phases imported with "
                  f"expected space group and cell.", flush=True)

        # ── Diagnostic: HAP keys (per-phase histogram-and-phase data) ─────
        # Used to verify the structure GSAS-II expects when we set per-phase
        # Size / Mustrain refinement flags.  Printed once at startup so
        # we can confirm the dict layout matches the installed GSAS-II.
        print("  HAP keys diagnostic (one-time, for per-phase refinement):",
              flush=True)
        for _pi, _po in enumerate(gsas_phases):
            try:
                _hap_dict = _po.data.get('Histograms', {})
                if not _hap_dict:
                    print(f"    Phase {_pi} ({_po.name}): no HAP entries.",
                          flush=True)
                    continue
                _hap = list(_hap_dict.values())[0]
                _keys = sorted(_hap.keys()) if isinstance(_hap, dict) else []
                print(f"    Phase {_pi} ({_po.name}): HAP keys = {_keys}",
                      flush=True)
                if isinstance(_hap, dict):
                    if 'Size' in _hap:
                        print(f"      Size:     {_hap['Size']}", flush=True)
                    if 'Mustrain' in _hap:
                        print(f"      Mustrain: {_hap['Mustrain']}",
                              flush=True)
                    if 'Pref.Ori.' in _hap:
                        print(f"      Pref.Ori.: {_hap['Pref.Ori.']}",
                              flush=True)
            except Exception as _e:
                print(f"    Phase {_pi}: HAP inspection failed: {_e}",
                      flush=True)

        # Track which stage succeeded (for fallback on failure)
        last_good_stage = 0
        _failed_stages = []

        # ── Helper to safely run a refinement stage ──────────────────────
        def _safe_refine(stage_name, refinement_dicts, stage_num):
            nonlocal last_good_stage
            try:
                gpx.do_refinements(refinement_dicts)
                last_good_stage = stage_num
                # Guard: ensure histogram scale stays fixed at 1.0.
                # GSAS-II's do_refinements can re-enable it internally
                # when 'Scale' appears in any refinement dict.
                try:
                    histogram.data['Sample Parameters']['Scale'] = [1.0, False]
                except (KeyError, TypeError):
                    pass
                return True
            except Exception as e:
                # GSAS-II internally restores from .bak0.gpx on failure,
                # so the project state reverts to pre-refinement. Safe to continue.
                _failed_stages.append(stage_name)
                warnings.warn(f"GSAS-II: {stage_name} failed ({e}). "
                             f"Continuing with results from stage {last_good_stage}.")
                return False

        # Determine structural complexity → boost cycles for complex phases
        def _count_asym_sites(ph):
            result = _reduce_to_asymmetric_unit(
                ph.get('cif_text', ''),
                declared_sg=ph.get('spacegroup_number'))
            if isinstance(result, str):
                # CifWriter returned full CIF text — parse site count
                _parsed = parse_cif(result)
                return len(_parsed.get('sites') or [])
            return len(result) if result else 0
        max_asym_atoms = max(
            (_count_asym_sites(ph) for ph in phases if ph.get('cif_text')),
            default=1)
        _complex = max_asym_atoms > 6  # W2C has ~6 asym sites
        _cyc_mult = 2 if _complex else 1  # double cycles for complex phases

        # ── Background coefficient selection ─────────────────────────────
        # Priority order:
        #   1. background_mode='fixedN' from options → exactly N terms
        #   2. auto_bg=False (user picked a number in GUI) → use n_bg_coeffs as-is
        #   3. auto_bg=True AND background_mode='auto' → DW/RMS auto-select
        #
        # The key distinction: when the user explicitly selects "6" or "8"
        # in the GUI dropdown, __init__.py sets auto_bg=False and n_bg_coeffs
        # to that value.  We must NOT override it with auto-selection.
        if background_mode.startswith('fixed'):
            try:
                n_bg_coeffs = int(background_mode.replace('fixed', ''))
                print(f"  Background: fixed mode → {n_bg_coeffs} Chebyshev "
                      f"coefficients.", flush=True)
            except ValueError:
                n_bg_coeffs = 6
        elif not auto_bg:
            # User explicitly selected a number in the GUI — respect it.
            print(f"  Background: user-specified {n_bg_coeffs} Chebyshev "
                  f"coefficients.", flush=True)
        elif auto_bg:
            n_bg_coeffs = _auto_select_bg_coeffs(
                tt, y_obs, phases, wavelength, tt_min, tt_max,
                user_n=n_bg_coeffs)
            print(f"  Background: auto-selected {n_bg_coeffs} Chebyshev "
                  f"coefficients.", flush=True)

        # Update the histogram's background config to match n_bg_coeffs.
        # This runs for ALL paths (auto, user-specified, fixed) so the
        # coefficient list is always the right length.
        bkg_data = histogram.data['Background']
        existing = bkg_data[0]
        old_n = existing[2] if len(existing) > 2 else 6
        existing[2] = n_bg_coeffs
        n_existing_coeffs = len(existing) - 3
        if n_bg_coeffs > n_existing_coeffs:
            existing.extend([0.0] * (n_bg_coeffs - n_existing_coeffs))
        elif n_bg_coeffs < n_existing_coeffs:
            bkg_data[0] = existing[:3 + n_bg_coeffs]

        if progress_callback:
            progress_callback('GSAS-II: stage 1 — refining background + scale...')

        # ── Stage 1: Background + scale ────────────────────────────────
        # IMPORTANT: Use EQUAL initial scales for all phases to avoid biasing
        # the refinement toward whichever phase is refined first.
        # Previous sequential approach (one phase at a time) caused the first
        # phase (WC) to absorb all intensity, driving the second phase (W2C)
        # scale to zero — making all subsequent stages fail for W2C.
        peak_height = float(np.max(y_r) - np.percentile(y_r, 5))
        n_phases = len(gsas_phases)
        init_scale = max(0.01, peak_height / max(1.0, n_phases * 50.0))
        print(f"  Initial scale for all {n_phases} phases: {init_scale:.4f} "
              f"(peak_height={peak_height:.1f})", flush=True)
        for phase_obj in gsas_phases:
            hapData = list(phase_obj.data['Histograms'].values())[0]
            hapData['Scale'] = [init_scale, False]  # start fixed

        # Stage 1a: refine phase scales FIRST with background fixed.
        # Prevents the background polynomial from absorbing real peak
        # intensity before the scale factors have a chance to activate.
        # This is especially important for highly crystalline patterns
        # (WC/W2C) where intense peaks sit on a relatively flat baseline.
        for phase_obj in gsas_phases:
            hapData = list(phase_obj.data['Histograms'].values())[0]
            hapData['Scale'] = [hapData['Scale'][0], True]  # turn on
        _safe_refine('scales only (background fixed)', [{
            'set': {},
            'cycles': min(max_cycles, 3),
        }], 1)

        # Stage 1b: refine background + scales together.
        # Now that scales are in the right ballpark, the background
        # polynomial can adjust without absorbing peak intensity.
        _safe_refine('background + scales', [{
            'set': {
                'Background': {'type': 'chebyschev-1', 'refine': True,
                                'no. coeffs': n_bg_coeffs},
            },
            'cycles': min(max_cycles, 5),
        }], 1)

        # ── Scale floor protection ──────────────────────────────────────
        # If any phase's scale collapsed to near-zero, it means the model
        # couldn't distribute intensity properly (perhaps wrong starting
        # guess).  Reset the collapsed phase to a floor value and re-refine
        # with background to give it another chance.
        _max_scale = max(
            list(po.data['Histograms'].values())[0].get('Scale', [0.001])[0]
            for po in gsas_phases)
        _scale_floor = max(1e-3, _max_scale * 0.01)  # at least 1% of max
        _any_reset = False
        for idx, phase_obj in enumerate(gsas_phases):
            hapData = list(phase_obj.data['Histograms'].values())[0]
            current_scale = hapData['Scale'][0]
            print(f"  Phase {idx} scale after simultaneous fit: "
                  f"{current_scale:.4e}", flush=True)
            if current_scale < _scale_floor:
                warnings.warn(
                    f"GSAS-II: phase {idx} scale ({current_scale:.2e}) below "
                    f"floor ({_scale_floor:.2e}) — resetting to floor value "
                    f"and re-refining.")
                hapData['Scale'] = [_scale_floor, True]
                _any_reset = True
        if _any_reset:
            # Re-refine with corrected scales + background
            _safe_refine('scales (after floor reset)', [{
                'set': {
                    'Background': {'type': 'chebyschev-1', 'refine': True,
                                    'no. coeffs': n_bg_coeffs},
                },
                'cycles': min(max_cycles, 8),
            }], 1)
            # Print scales after re-refinement
            for idx, phase_obj in enumerate(gsas_phases):
                hapData = list(phase_obj.data['Histograms'].values())[0]
                print(f"  Phase {idx} scale after re-refinement: "
                      f"{hapData['Scale'][0]:.4e}", flush=True)

        if progress_callback:
            progress_callback('GSAS-II: stage 2 — refining profile parameters...')

        # ── Stage 2: Profile + background + zero + transparency ──────
        # Refine Gaussian (U, W) and Lorentzian (X, Y) profile parameters
        # together with zero shift, background, and sample displacement.
        #
        # Profile parameter selection for lab XRD stability:
        #   Refined: U, W (Gaussian Caglioti), X, Y (Lorentzian), Zero,
        #            Displacement (geometry-dependent: DisplaceY for
        #            capillary/transmission, DisplaceX for Bragg-Brentano)
        #   Fixed:   V (99%+ correlated with U and W — instrument-determined)
        #            SH/L (instrument constant, 99.9% correlated with Zero)
        _displace_stg2 = ('DisplaceX' if geometry == 'bragg_brentano'
                          else 'DisplaceY')
        # When a measured .instprm is loaded, U/V/W/X/Y are FIXED —
        # only Zero and displacement refine alongside background.
        # This is the whole point of the measured standard: instrument
        # broadening is pinned, leaving only sample broadening for
        # Size/Mustrain to capture.
        if _measured_instprm:
            # With measured instprm: U/V/W (Gaussian Caglioti) and X
            # (Lorentzian strain) are fixed at the Si-standard values —
            # they're instrument-determined.  Zero is also fixed; all
            # position correction goes through DisplaceX/Y.
            #
            # Y (Lorentzian, FWHM ~ Y/cos θ) STAYS REFINABLE because it
            # absorbs sample size broadening.  A Si standard pins the
            # instrumental Y; nanocrystalline samples need additional
            # Y on top of that.  Locking Y at the Si value would force
            # nano peaks to look like bulk-Si peaks → broken fit for
            # WC/W2C (~9 nm crystallites).
            _inst_params_stg2 = ['Y']
            _stg2_label = 'bg + displacement + Y (U/V/W/X/SH/L/Zero fixed by instprm)'
            print(f"  Stage 2: U/V/W/X/SH/L/Zero fixed at instprm values; "
                  f"Y refinable for sample size broadening.", flush=True)
        elif verification_mode:
            # No measured instprm + verification_mode: tighten the profile
            # to reduce X↔Y correlation.  X is fixed at 0 (strain second-
            # order for nanocrystals).  Zero is fixed at the seed value;
            # DisplaceX/Y absorbs position correction.  V is fixed
            # (correlated with U and W).  This leaves U, W, Y refinable —
            # the minimum set that can describe Caglioti+Lorentzian width.
            _inst_params_stg2 = ['U', 'W', 'Y']
            _stg2_label = ('verify: bg + displacement + U + W + Y '
                           '(V/X/Zero/SH/L fixed)')
            print(f"  Stage 2 (verification_mode, no instprm): "
                  f"refining U, W, Y; V/X/Zero/SH/L held fixed.",
                  flush=True)
        else:
            _inst_params_stg2 = ['U', 'W', 'X', 'Y', 'Zero']
            _stg2_label = 'profile + bg + zero + displacement'

        # Swap Zero ↔ Displace (when verify_use_zero_not_displace is on):
        # Add Zero to the refinable list, drop displacement from the
        # Sample Parameters list, and explicitly zero the DisplaceX/Y
        # value so it can't drift from a previous stage.
        if use_zero_not_displace_opt:
            if 'Zero' not in _inst_params_stg2:
                _inst_params_stg2 = list(_inst_params_stg2) + ['Zero']
            try:
                for _disp_key in ('DisplaceX', 'DisplaceY'):
                    _dp = histogram.data['Sample Parameters'].get(_disp_key)
                    if _dp and isinstance(_dp, list):
                        _dp[0] = 0.0
                        if len(_dp) >= 2:
                            _dp[1] = False
            except Exception as _e:
                print(f"    Warning: could not pin displacement to 0: {_e}",
                      flush=True)
            _stg2_sample_params = []
            _stg2_label = (_stg2_label.rstrip('.')
                           + ' [Zero ON, Displace=0 fixed]')
        else:
            _stg2_sample_params = [_displace_stg2]

        # Diagnostic: add X (Lorentzian strain) to the refine list.
        # Stacks on top of any of the above branches.
        if refine_x_opt and 'X' not in _inst_params_stg2:
            _inst_params_stg2 = list(_inst_params_stg2) + ['X']
            _stg2_label = _stg2_label + ' +X'

        # Fix-Y dominates: drop Y from refine list and pin its value.
        if fix_y_opt:
            _inst_params_stg2 = [p for p in _inst_params_stg2 if p != 'Y']
            try:
                _inst_dict = histogram.data['Instrument Parameters'][0]
                if 'Y' in _inst_dict and isinstance(_inst_dict['Y'], list):
                    _y_target = (float(y_fixed_value_opt)
                                 if y_fixed_value_opt is not None
                                 else float(_inst_dict['Y'][1]))
                    _inst_dict['Y'][1] = _y_target  # current value
                    if len(_inst_dict['Y']) >= 3:
                        _inst_dict['Y'][2] = False   # refine flag off
                    print(f"    Stage 2: Y fixed at {_y_target:.4f}.",
                          flush=True)
            except Exception as _e:
                print(f"    Warning: could not fix Y: {_e}", flush=True)
            _stg2_label = _stg2_label + ' [Y fixed]'

        _safe_refine(_stg2_label, [{
            'set': {
                'Background': {'type': 'chebyschev-1', 'refine': True,
                                'no. coeffs': n_bg_coeffs},
                'Instrument Parameters': _inst_params_stg2,
                'Sample Parameters': _stg2_sample_params,
            },
            'cycles': min(max_cycles, 10 * _cyc_mult),
        }], 2)

        # ── Identify phases with negligible scale ──────────────────────────
        # Phases with effectively zero scale contribute nothing to the pattern.
        # Refining their cell/atom parameters has zero gradient and causes
        # divergence (e.g. "Invalid cell metric tensor" errors).  Skip them.
        _negligible_phases = set()
        _max_scale_now = max(
            list(po.data['Histograms'].values())[0].get('Scale', [0])[0]
            for po in gsas_phases) or 1e-10
        for idx, phase_obj in enumerate(gsas_phases):
            hapData = list(phase_obj.data['Histograms'].values())[0]
            s = hapData.get('Scale', [0])[0]
            if s < _max_scale_now * 1e-3:  # less than 0.1% of max
                _negligible_phases.add(idx)
                print(f"  Phase {idx}: scale={s:.2e} << max={_max_scale_now:.2e} "
                      f"— skipping cell/atom refinement stages.", flush=True)

        # ── Preferred orientation setup ──────────────────────────────────
        # March-Dollase correction for preferred orientation.  Controlled
        # by options['preferred_orientation']:
        #   "auto"  → enable for hexagonal/trigonal phases only (default)
        #   "off"   → disable for all phases
        #   "force" → enable for all phases
        #
        # Common preferred orientation directions:
        #   Hexagonal platelets (WC, BN, MoS2):  [0, 0, 1]
        #   Orthorhombic needles:                 [1, 0, 0] or [0, 1, 0]
        #
        # The MD ratio starts at 1.0 (no texture) and is refined later
        # alongside Uiso (both affect peak intensities).
        _pref_ori_phases = set()
        # Track which phases want PO held at a fixed value (per-phase).
        # Used by the refine-flag loop later to skip setting Pref.Ori.→True.
        _po_fixed_phases = set()
        _po_phase_modes = {}
        for idx, (phase_obj, ph_input) in enumerate(zip(gsas_phases, phases)):
            if idx in _negligible_phases:
                continue
            sys_ = (ph_input.get('system') or 'triclinic').lower()
            _ph_name = getattr(phase_obj, 'name', '') or ''
            _popts = _phase_opts_for(idx, _ph_name)

            # Per-phase override path: phase_options entry wins.
            _po_mode_p = _popts.get('po_mode', 'off')
            _po_axis_p = _popts.get('po_axis')
            _po_value_p = _popts.get('po_value', 1.0)
            _has_phase_options_entry = (
                idx < len(phase_options_list)
                and isinstance(phase_options_list[idx], dict))
            _has_per_phase_po = (_po_mode_p in ('fixed', 'refined'))

            # Legacy global path: preferred_orientation + fix_po_opt.
            _enable_md_global = False
            if (not _has_phase_options_entry and not _has_per_phase_po
                    and preferred_orientation != 'off'):
                if preferred_orientation == 'force':
                    _enable_md_global = True
                elif preferred_orientation == 'auto':
                    _enable_md_global = sys_ in ('hexagonal', 'trigonal')

            if not (_has_per_phase_po or _enable_md_global):
                continue

            # Resolve direction: explicit per-phase axis wins, else
            # crystal-system default.
            if _po_axis_p and isinstance(_po_axis_p, (list, tuple)) \
                    and len(_po_axis_p) == 3:
                _po_dir = [int(round(float(x))) for x in _po_axis_p]
            else:
                if sys_ in ('hexagonal', 'trigonal'):
                    _po_dir = [0, 0, 1]
                elif sys_ == 'orthorhombic':
                    _po_dir = [1, 0, 0]
                else:
                    _po_dir = [0, 0, 1]

            # Resolve value: per-phase value wins, else fix_po_opt value
            # (legacy), else 1.0.
            if _has_per_phase_po:
                _po_init = float(_po_value_p) if _po_mode_p == 'fixed' else 1.0
                _is_fixed = (_po_mode_p == 'fixed')
            else:
                _po_init = (po_fixed_value_opt if fix_po_opt else 1.0)
                _is_fixed = fix_po_opt

            try:
                hapData = list(phase_obj.data['Histograms'].values())[0]
                po = hapData.get('Pref.Ori.', ['MD', 1.0, False, _po_dir])
                po[0] = 'MD'
                po[1] = _po_init
                po[3] = _po_dir
                hapData['Pref.Ori.'] = po
                _pref_ori_phases.add(idx)
                _po_phase_modes[idx] = (
                    _po_mode_p if _has_per_phase_po else 'global')
                if _is_fixed:
                    _po_fixed_phases.add(idx)
                _src = 'phase_options' if _has_per_phase_po else 'auto/global'
                if _is_fixed:
                    print(f"  Phase {idx} ({sys_}): March-Dollase "
                          f"{_po_dir} PO held at "
                          f"{_po_init:.4f} (not refined) [{_src}].",
                          flush=True)
                else:
                    print(f"  Phase {idx} ({sys_}): March-Dollase "
                          f"{_po_dir} preferred orientation enabled "
                          f"[{_src}].", flush=True)
            except Exception as e:
                print(f"  Phase {idx}: could not set preferred "
                      f"orientation: {e}", flush=True)
        if not _pref_ori_phases:
            print("  Preferred orientation: no phases enabled.", flush=True)

        if progress_callback:
            progress_callback('GSAS-II: stage 3 — refining cell parameters...')

        # ── Stage 3: Cell parameters (one phase at a time) ───────────────
        # Refining all cells at once with many atoms can cause arccos
        # errors when cell angles go unphysical. Do it per phase.
        # For cubic phases, lock angles to 90° to prevent arccos errors.
        if verification_mode:
            print("  Stage 3 (cell): skipped (verification_mode).", flush=True)
        for idx, (phase_obj, ph_input) in enumerate(zip(gsas_phases, phases)):
            if verification_mode:
                continue
            if idx in _negligible_phases:
                print(f"  Skipping cell refinement for phase {idx} "
                      f"(negligible scale).", flush=True)
                continue
            sys_ = (ph_input.get('system') or 'triclinic').lower()
            # For high-symmetry systems, clamp angles before refining
            # to prevent them from drifting to unphysical values.
            if sys_ in ('cubic', 'tetragonal', 'orthorhombic'):
                try:
                    _cell = phase_obj.data['General']['Cell']
                    _cell[4] = 90.0  # alpha
                    _cell[5] = 90.0  # beta
                    _cell[6] = 90.0  # gamma
                except Exception:
                    pass
            elif sys_ in ('hexagonal', 'trigonal'):
                try:
                    _cell = phase_obj.data['General']['Cell']
                    _cell[4] = 90.0   # alpha
                    _cell[5] = 90.0   # beta
                    _cell[6] = 120.0  # gamma
                except Exception:
                    pass
            # Save cell params before refinement for sanity check
            _pre_cell = list(phase_obj.data['General']['Cell'])
            phase_obj.set_refinements({'Cell': True})
            ok = _safe_refine(f'cell (phase {idx})', [{
                'set': {},
                'cycles': min(max_cycles, 8 * _cyc_mult),
            }], 3)
            if not ok:
                # Turn cell refinement back off for this phase
                phase_obj.set_refinements({'Cell': False})
            else:
                # ── Post-stage cell sanity check ──────────────────────
                # Even when _safe_refine succeeds, the cell may have
                # diverged to unphysical values within the cycles before
                # GSAS-II's internal checks caught it. Detect and revert.
                _post_cell = phase_obj.data['General']['Cell']
                _cell_ok = True
                for _ci in range(1, 4):  # a, b, c
                    _pre_v = float(_pre_cell[_ci])
                    _post_v = float(_post_cell[_ci])
                    if _pre_v > 0.1:
                        _ratio = _post_v / _pre_v
                        if _ratio > 2.0 or _ratio < 0.5:
                            _cell_ok = False
                            break
                for _ci in range(4, 7):  # alpha, beta, gamma
                    if abs(float(_post_cell[_ci]) - float(_pre_cell[_ci])) > 30:
                        _cell_ok = False
                        break
                if not _cell_ok:
                    warnings.warn(
                        f"GSAS-II: cell params for phase {idx} diverged "
                        f"after refinement — reverting to pre-refinement "
                        f"values and turning off cell refinement.")
                    for _ci in range(len(_pre_cell)):
                        phase_obj.data['General']['Cell'][_ci] = _pre_cell[_ci]
                    phase_obj.set_refinements({'Cell': False})

        if progress_callback:
            progress_callback('GSAS-II: stage 4 — refining atomic displacement + size/strain...')

        # ── Stage 4: Atomic displacement (Uiso) + preferred orientation ──
        #             + Size and Mustrain HAP parameters
        # Uiso and preferred orientation both affect peak intensities, so
        # they should refine together to avoid one absorbing the other's
        # effect.
        #
        # NEW: Also refine isotropic Size and Mustrain (crystallite size
        # and micro-strain) per-phase.  These HAP parameters control
        # per-phase Lorentzian broadening in GSAS-II's profile model.
        # Without refining them, the instrument Y parameter (shared across
        # all phases) must accommodate all size broadening, which fails
        # for mixed-phase systems (e.g. WC + W2C) where phases have
        # different crystallite sizes.  Unrefined Size defaults to
        # GSAS-II's initial value (~1 µm), producing wrong peak shapes
        # and trapping the refinement in a local minimum.
        # Stage 4 runs when not in verification_mode, OR when explicitly
        # opted into Uiso refinement via verify_refine_uiso.
        _run_stage4 = (not verification_mode) or refine_uiso_opt
        if not _run_stage4:
            print("  Stage 4 (Uiso + MD): skipped (verification_mode, "
                  "no verify_refine_uiso).", flush=True)
        elif verification_mode and refine_uiso_opt:
            print("  Stage 4 (Uiso): running in verification_mode "
                  "(verify_refine_uiso=True).", flush=True)
        for idx, phase_obj in enumerate(gsas_phases):
            if not _run_stage4:
                continue
            if idx in _negligible_phases:
                continue
            try:
                phase_obj.set_refinements({'Atoms': {'all': 'U'}})
            except Exception:
                pass
        # Enable preferred orientation refinement for applicable phases.
        # Per-phase: phases whose PO is "fixed" (in phase_options or via
        # legacy fix_po_opt) keep refine_flag=False.  Others get the
        # refine flag set when verification mode allows it.
        for idx in _pref_ori_phases:
            _po_mode_for_idx = _po_phase_modes.get(idx, 'global')
            if (verification_mode and not verify_refine_po_opt
                    and _po_mode_for_idx != 'refined'):
                continue
            if idx in _po_fixed_phases:
                # Defensive: keep refine flag cleared.
                if idx not in _negligible_phases:
                    try:
                        phase_obj_p = gsas_phases[idx]
                        _hap_p = list(
                            phase_obj_p.data['Histograms'].values())[0]
                        _po_p = _hap_p.get('Pref.Ori.')
                        if _po_p and len(_po_p) >= 3:
                            _po_p[2] = False
                    except Exception:
                        pass
                continue
            if idx not in _negligible_phases:
                try:
                    phase_obj = gsas_phases[idx]
                    phase_obj.set_HAP_refinements(
                        {'Pref.Ori.': True})
                    _po_via = ('Stage 6' if (verification_mode
                                              and not _run_stage4)
                               else 'Stage 4')
                    print(f"  Phase {idx}: March-Dollase ratio refinement "
                          f"enabled (refines in {_po_via}).", flush=True)
                except Exception as e:
                    print(f"  Phase {idx}: could not enable Pref.Ori. "
                          f"refinement: {e}", flush=True)
        if _run_stage4:
            _safe_refine('Uiso + preferred orientation', [{
                'set': {},
                'cycles': min(max_cycles, 8 * _cyc_mult),
            }], 4)

        # ── Stage 4b: Crystallite Size + Microstrain (per-phase) ────────
        # Drive per-phase HAP refinement from phase_options (or the
        # legacy WC/W2C name-based fallback inside _phase_opts_for).
        # No phase-name assumptions in this loop — it just iterates
        # phase indices and applies whatever options resolve.

        # Pre-compute resolved options per phase
        _resolved_phase_opts = []
        for idx, phase_obj in enumerate(gsas_phases):
            _resolved_phase_opts.append(
                _phase_opts_for(idx, getattr(phase_obj, 'name', '')))

        _any_per_phase_hap = any(
            (po.get('refine_size') or po.get('refine_mustrain'))
            for po in _resolved_phase_opts)
        _run_stage4b = ((not verification_mode) or refine_size_opt
                        or _any_per_phase_hap)
        if not _run_stage4b:
            print("  Stage 4b (HAP size/mustrain): skipped "
                  "(verification_mode, no per-phase HAP options).",
                  flush=True)
        elif verification_mode:
            print("  Stage 4b (HAP size/mustrain): running in "
                  "verification_mode.", flush=True)

        # Per-phase decision and refinement-flag setting
        _stage4b_anything_set = False
        for idx, phase_obj in enumerate(gsas_phases):
            if not _run_stage4b:
                continue
            if idx in _negligible_phases:
                continue
            _popts = _resolved_phase_opts[idx]
            _enable_size = bool(refine_size_opt or _popts.get('refine_size'))
            _enable_mustrain = bool(_popts.get('refine_mustrain'))
            _ph_label = getattr(phase_obj, 'name', f'phase{idx}')

            if _enable_size:
                try:
                    phase_obj.set_HAP_refinements(
                        {'Size': {'type': 'isotropic', 'refine': True}})
                    print(f"  Phase {idx} ({_ph_label}): Size refinement "
                          f"enabled (isotropic).", flush=True)
                    _stage4b_anything_set = True
                except Exception as e:
                    print(f"  Phase {idx}: could not enable Size "
                          f"refinement: {e}", flush=True)
            if _enable_mustrain:
                try:
                    phase_obj.set_HAP_refinements(
                        {'Mustrain': {'type': 'isotropic', 'refine': True}})
                    print(f"  Phase {idx} ({_ph_label}): Mustrain "
                          f"refinement enabled (isotropic).", flush=True)
                    _stage4b_anything_set = True
                except Exception as e:
                    print(f"  Phase {idx}: could not enable Mustrain "
                          f"refinement: {e}", flush=True)

        if _run_stage4b and _stage4b_anything_set:
            _safe_refine('HAP per-phase: size and/or mustrain', [{
                'set': {},
                'cycles': min(max_cycles, 8 * _cyc_mult),
            }], 4)
        elif _run_stage4b:
            print("  Stage 4b: no phases matched any HAP option, "
                  "skipping refinement call.", flush=True)

        # Check for instrument-limited Size (≥ 0.9 µm ≈ 9000 Å).
        # When Stage 4b was skipped, values are GSAS-II defaults and
        # meaningless, so report "not refined" instead.
        for idx, phase_obj in enumerate(gsas_phases):
            if idx in _negligible_phases:
                continue
            if not _run_stage4b:
                print(f"  Phase {idx}: Size: not refined "
                      f"(Stage 4b skipped).", flush=True)
                continue
            try:
                hapData = list(phase_obj.data['Histograms'].values())[0]
                size_data = hapData.get('Size', [])
                if size_data and len(size_data) > 1:
                    size_um = float(size_data[1][0])
                    if size_um >= 0.9:
                        print(f"  Phase {idx}: Size = {size_um:.3f} µm "
                              f"(instrument-limited, ≥ 1 µm).", flush=True)
                    else:
                        print(f"  Phase {idx}: Size = {size_um:.4f} µm "
                              f"= {size_um * 10000:.0f} Å.", flush=True)
            except Exception:
                pass

        # Freeze Uiso after Stage 4.  For heavy elements like W (Z=74)
        # with limited lab XRD reflections, Uiso is weakly determined
        # and becomes singular when refined alongside scale + profile in
        # later stages.  The Stage 4 value is our best estimate; freezing
        # it prevents SVD drops in Stages 5–6.
        for phase_obj in gsas_phases:
            try:
                phase_obj.clear_refinements({'Atoms': {'all': ''}})
            except Exception:
                pass

        # ── Stage 5: Atom positions (XYZ) ─────────────────────────────
        # Gated by options['refine_xyz'] (default False).
        # Without a measured .instprm, XYZ refinement can absorb errors
        # from polarization, absorption, PO, background, or wrong phase
        # choice — making Rwp look better while hiding real problems.
        # Enable only when instrument params are pinned.
        if refine_xyz:
            if progress_callback:
                progress_callback('GSAS-II: stage 5 — refining atom positions (XYZ)...')

            # Save atom positions before refinement for damping check.
            saved_xyz = {}
            for idx, phase_obj in enumerate(gsas_phases):
                if idx in _negligible_phases:
                    continue  # skip XYZ for phases with negligible scale
                phase_atoms = phase_obj.data.get('Atoms', [])
                saved_xyz[idx] = [(a[3], a[4], a[5]) if len(a) > 5 else None
                                  for a in phase_atoms]
                try:
                    phase_obj.set_refinements({'Atoms': {'all': 'X'}})
                except Exception:
                    pass
            xyz_ok = _safe_refine('atom XYZ', [{
                'set': {},
                'cycles': min(max_cycles, 3 * _cyc_mult),
            }], 5)

            # Damping check: if any atom jumped > 0.5 fractional units,
            # it likely hopped to a symmetry-equivalent site — revert it.
            if xyz_ok:
                for idx, phase_obj in enumerate(gsas_phases):
                    phase_atoms = phase_obj.data.get('Atoms', [])
                    for j, atom in enumerate(phase_atoms):
                        if (len(atom) > 5 and idx in saved_xyz
                                and j < len(saved_xyz[idx])
                                and saved_xyz[idx][j] is not None):
                            ox, oy, oz = saved_xyz[idx][j]
                            dx = abs(atom[3] - ox)
                            dy = abs(atom[4] - oy)
                            dz = abs(atom[5] - oz)
                            if max(dx, dy, dz) > 0.5:
                                warnings.warn(
                                    f"GSAS-II: atom {j} in phase {idx} jumped "
                                    f"by ({dx:.3f},{dy:.3f},{dz:.3f}) — "
                                    f"reverting to original position.")
                                atom[3], atom[4], atom[5] = ox, oy, oz

            # Turn off all atom refinement flags
            for phase_obj in gsas_phases:
                try:
                    phase_obj.clear_refinements({'Atoms': {'all': ''}})
                except Exception:
                    pass
        else:
            print("  Stage 5: XYZ refinement skipped (refine_xyz=False).",
                  flush=True)

        # Ensure atom refinement flags are off regardless of XYZ path
        for phase_obj in gsas_phases:
            try:
                phase_obj.clear_refinements({'Atoms': {'all': ''}})
            except Exception:
                pass

        if progress_callback:
            progress_callback('GSAS-II: stage 6 — final all-together refinement...')

        # ── Stage 6: Final all-together refinement ──────────────────────
        # Refine background, profile (U, W, X, Y, Zero), scales, and cell
        # parameters simultaneously.  Uiso is frozen at the Stage 4 value
        # to avoid SVD singularity (Uiso ↔ scale correlation is too strong
        # for heavy elements like W with limited lab XRD reflections).
        # Cell is re-enabled; Uiso stays frozen (determined in Stage 4).
        #
        # Cell-refinement gate for Stage 6:
        #   verification_mode=True, verify_refine_cell=False  → Cell fixed
        #   verification_mode=True, verify_refine_cell=True   → Cell on
        #   verification_mode=False                            → Cell on (default)
        if verification_mode and not verify_refine_cell:
            print("  Stage 6: cell refinement disabled (verification_mode).",
                  flush=True)
            for idx, phase_obj in enumerate(gsas_phases):
                try:
                    phase_obj.set_refinements({'Cell': False})
                except Exception:
                    pass
        else:
            if verification_mode and verify_refine_cell:
                print("  Stage 6: cell refinement enabled "
                      "(verify_refine_cell=True).", flush=True)
            for idx, phase_obj in enumerate(gsas_phases):
                if idx not in _negligible_phases:
                    phase_obj.set_refinements({'Cell': True})

        # ── Snapshot W2C cell BEFORE Stage 6 (for uniform-scale post-fit) ──
        _uniform_scale_snapshot = {}    # idx → (a0, b0, c0, V0)
        if uniform_cell_w2c_opt and verify_refine_cell:
            for _idx_s, (_phase_obj_s, _ph_input_s) in enumerate(
                    zip(gsas_phases, phases)):
                if _uniform_cell_phase_idxs and _idx_s not in _uniform_cell_phase_idxs:
                    continue
                _sg_n = int(_ph_input_s.get('spacegroup_number', 0) or 0)
                _formula = str(_ph_input_s.get('formula', '')).replace(' ', '')
                _is_pbcn_w2c = (_sg_n == 60
                                and ('W2C' in _formula or 'CW2' in _formula))
                if not _is_pbcn_w2c:
                    continue
                try:
                    _cell_arr = _phase_obj_s.data['General']['Cell']
                    _a0 = float(_cell_arr[1])
                    _b0 = float(_cell_arr[2])
                    _c0 = float(_cell_arr[3])
                    _V0 = _a0 * _b0 * _c0
                    _uniform_scale_snapshot[_idx_s] = (_a0, _b0, _c0, _V0)
                    print(f"  Uniform-scale W2C: snapshot phase {_idx_s} "
                          f"a={_a0:.5f}, b={_b0:.5f}, c={_c0:.5f}, "
                          f"V={_V0:.4f} Å³", flush=True)
                except Exception as _e:
                    print(f"  Uniform-scale W2C: snapshot failed for "
                          f"phase {_idx_s}: {_e}", flush=True)
        # Displacement correction depends on geometry:
        #   capillary       → DisplaceY (radial displacement from beam axis)
        #   bragg_brentano  → DisplaceX (sample-height/flat-plate displacement)
        # Using the wrong model distorts peak positions and lattice params.
        _displace_param = 'DisplaceX' if geometry == 'bragg_brentano' else 'DisplaceY'
        print(f"  Stage 6: displacement model = {_displace_param} "
              f"(geometry={geometry})", flush=True)

        # With measured .instprm: fix U/V/W/X (instrument-determined) and
        # Zero/SH/L (geometric).  Y stays refinable to absorb sample size
        # broadening — see Stage 2 comment.
        if _measured_instprm:
            _inst_params_stg6 = ['Y']
            _stg6_label = ('final: cell + bg + displacement + Y '
                           '(U/V/W/X/SH/L/Zero fixed by instprm)')
        elif verification_mode:
            # No measured instprm + verification_mode: same restricted
            # profile as Stage 2 — U, W, Y only.  Cell is gated separately
            # by verify_refine_cell (cell-flag block above).
            _inst_params_stg6 = ['U', 'W', 'Y']
            _cell_state = 'cell on' if verify_refine_cell else 'cell fixed'
            _stg6_label = (f'verify-final: bg + displacement + U + W + Y '
                           f'({_cell_state}; V/X/Zero/SH/L fixed)')
        else:
            _inst_params_stg6 = ['U', 'W', 'X', 'Y', 'Zero']
            _stg6_label = 'final: all parameters together'

        # Swap Zero ↔ Displace for Stage 6 (same as Stage 2)
        if use_zero_not_displace_opt:
            if 'Zero' not in _inst_params_stg6:
                _inst_params_stg6 = list(_inst_params_stg6) + ['Zero']
            try:
                for _disp_key in ('DisplaceX', 'DisplaceY'):
                    _dp = histogram.data['Sample Parameters'].get(_disp_key)
                    if _dp and isinstance(_dp, list):
                        _dp[0] = 0.0
                        if len(_dp) >= 2:
                            _dp[1] = False
            except Exception as _e:
                print(f"    Warning: could not pin displacement to 0 "
                      f"in Stage 6: {_e}", flush=True)
            _stg6_sample_params = []
            _stg6_label = (_stg6_label.rstrip(')').rstrip('.')
                           + ' [Zero ON, Displace=0 fixed]')
        else:
            _stg6_sample_params = [_displace_param]

        # Diagnostic: add X (Lorentzian strain) to the refine list.
        if refine_x_opt and 'X' not in _inst_params_stg6:
            _inst_params_stg6 = list(_inst_params_stg6) + ['X']
            _stg6_label = _stg6_label + ' +X'

        # Fix-Y for Stage 6 (matches Stage 2 logic)
        if fix_y_opt:
            _inst_params_stg6 = [p for p in _inst_params_stg6 if p != 'Y']
            try:
                _inst_dict_6 = histogram.data['Instrument Parameters'][0]
                if 'Y' in _inst_dict_6 and isinstance(_inst_dict_6['Y'], list):
                    _y_target_6 = (float(y_fixed_value_opt)
                                    if y_fixed_value_opt is not None
                                    else float(_inst_dict_6['Y'][1]))
                    _inst_dict_6['Y'][1] = _y_target_6
                    if len(_inst_dict_6['Y']) >= 3:
                        _inst_dict_6['Y'][2] = False
                    print(f"    Stage 6: Y fixed at {_y_target_6:.4f}.",
                          flush=True)
            except Exception as _e:
                print(f"    Warning: could not fix Y in Stage 6: {_e}",
                      flush=True)
            _stg6_label = _stg6_label + ' [Y fixed]'

        _safe_refine(_stg6_label, [{
            'set': {
                'Background': {'type': 'chebyschev-1', 'refine': True,
                                'no. coeffs': n_bg_coeffs},
                'Instrument Parameters': _inst_params_stg6,
                'Sample Parameters': _stg6_sample_params,
            },
            'cycles': min(max_cycles, 15 * _cyc_mult),
        }], 6)

        # ── Y non-negative enforcement ─────────────────────────────────
        # Negative Y is unphysical (Lorentzian width can't be negative)
        # but the LM optimizer can drive Y < 0 when X and Y trade off.
        # When verify_y_nonnegative is on, clamp refined Y to 0 if it
        # went negative, then run a 0-cycle refresh so ycalc updates.
        if y_nonneg_opt and not fix_y_opt:
            try:
                _inst_dict_y = histogram.data['Instrument Parameters'][0]
                if 'Y' in _inst_dict_y and isinstance(_inst_dict_y['Y'], list):
                    _y_now = float(_inst_dict_y['Y'][1])
                    if _y_now < 0.0:
                        print(f"  Y nonneg: refined Y = {_y_now:.4f} < 0; "
                              f"clamping to 0 and refreshing ycalc.",
                              flush=True)
                        _y_nonnegative_clamped = True
                        _inst_dict_y['Y'][1] = 0.0
                        if len(_inst_dict_y['Y']) >= 3:
                            _inst_dict_y['Y'][2] = False
                        # Turn off cell so the refresh doesn't drift
                        for _phase_obj_y in gsas_phases:
                            try:
                                _phase_obj_y.set_refinements({'Cell': False})
                            except Exception:
                                pass
                        _safe_refine('Y nonneg: refresh ycalc with Y=0',
                                     [{'set': {}, 'cycles': 0}], 6)
                    else:
                        print(f"  Y nonneg: refined Y = {_y_now:.4f} ≥ 0; "
                              f"no clamping needed.", flush=True)
            except Exception as _e:
                print(f"  Y nonneg: failed to inspect/clamp Y: {_e}",
                      flush=True)

        # ── Branch B: uniform-cell scaling for W2C ─────────────────────
        # Stage 6 refined (a, b, c) for W2C independently.  Now post-
        # process: compute s = (V_refined / V_start)^(1/3) and apply
        # (a, b, c) ← (s·a₀, s·b₀, s·c₀).  This forces the cell change
        # to be a uniform contraction (preserves axis ratios from
        # mp-2034) while keeping the volume that the data preferred.
        # Run a 0-cycle refinement so GSAS-II re-evaluates ycalc with
        # the uniformly-scaled cell.  χ² will be slightly worse than
        # the unconstrained fit if the data wanted real anisotropy, or
        # essentially the same if the anisotropy was data-resolution
        # noise.
        if uniform_cell_w2c_opt and verify_refine_cell \
                and _uniform_scale_snapshot:
            _scaled_any = False
            for idx_h, (a0, b0, c0, V0) in _uniform_scale_snapshot.items():
                if idx_h >= len(gsas_phases):
                    continue
                _phase_obj_h = gsas_phases[idx_h]
                try:
                    _cell = _phase_obj_h.data['General']['Cell']
                    a_ref = float(_cell[1])
                    b_ref = float(_cell[2])
                    c_ref = float(_cell[3])
                    V_ref = a_ref * b_ref * c_ref
                    if V0 <= 0:
                        print(f"  Uniform-scale W2C: phase {idx_h} has "
                              f"non-positive starting volume; skipping.",
                              flush=True)
                        continue
                    s = (V_ref / V0) ** (1.0 / 3.0)
                    a_new = s * a0
                    b_new = s * b0
                    c_new = s * c0
                    print(f"  Uniform-scale W2C: phase {idx_h} — "
                          f"refined a={a_ref:.5f}, b={b_ref:.5f}, "
                          f"c={c_ref:.5f} (V={V_ref:.4f}); applying "
                          f"s={s:.5f} → a={a_new:.5f}, b={b_new:.5f}, "
                          f"c={c_new:.5f} (V={a_new*b_new*c_new:.4f})",
                          flush=True)
                    _cell[1] = a_new
                    _cell[2] = b_new
                    _cell[3] = c_new
                    _cell[7] = a_new * b_new * c_new
                    _scaled_any = True
                except Exception as _e:
                    print(f"  Uniform-scale W2C: failed for phase "
                          f"{idx_h}: {_e}", flush=True)
            if _scaled_any:
                # Turn off cell refinement everywhere and run 0 cycles
                # so ycalc re-computes with the uniformly-scaled W2C cell.
                for _phase_obj_h in gsas_phases:
                    try:
                        _phase_obj_h.set_refinements({'Cell': False})
                    except Exception:
                        pass
                _safe_refine('uniform-scale: refresh ycalc with '
                             'uniformly-scaled W2C cell',
                             [{'set': {}, 'cycles': 0}], 6)

        # ── Post-Stage-6 diagnostic: actual GSAS-II parameter state ────
        # The phase_results table reports cooked/rounded numbers from
        # extracted dicts.  Print the raw sample + instrument parameter
        # entries here so we can see GSAS-II's actual stored values
        # (including refine flags) without trusting the display layer.
        print("  Sample parameters after Stage 6:", flush=True)
        try:
            _sample = histogram.data.get('Sample Parameters', {})
            for _k in sorted(_sample.keys()):
                if any(tok in _k for tok in ('Displace', 'Shift', 'Zero')):
                    print(f"    {_k}: {_sample[_k]}", flush=True)
        except Exception as _e:
            print(f"    Could not inspect sample parameters: {_e}",
                  flush=True)
        print("  Instrument parameters after Stage 6:", flush=True)
        try:
            _inst_dump = histogram.data.get('Instrument Parameters',
                                             [{}])[0]
            for _k in ('Zero', 'X', 'Y', 'U', 'V', 'W', 'SH/L'):
                if _k in _inst_dump:
                    print(f"    {_k}: {_inst_dump[_k]}", flush=True)
        except Exception as _e:
            print(f"    Could not inspect instrument parameters: {_e}",
                  flush=True)
        # Per-phase HAP scales and Pref.Ori. — confirms whether MD
        # ratio was actually written/refined.
        print("  Per-phase HAP state after Stage 6:", flush=True)
        try:
            for _idx_d, _phase_obj_d in enumerate(gsas_phases):
                _hapData = list(_phase_obj_d.data['Histograms'].values())[0]
                _scale = _hapData.get('Scale', [None])
                _po = _hapData.get('Pref.Ori.')
                _name = getattr(_phase_obj_d, 'name', f'Phase {_idx_d}')
                print(f"    Phase {_idx_d} ({_name}): "
                      f"Scale = {_scale}, Pref.Ori. = {_po}", flush=True)
        except Exception as _e:
            print(f"    Could not inspect HAP state: {_e}", flush=True)

        # Turn off cell refinement after final stage to avoid affecting
        # any subsequent operations
        for phase_obj in gsas_phases:
            try:
                phase_obj.set_refinements({'Cell': False})
            except Exception:
                pass

        # Warn if multiple refinement stages failed — result may be unreliable
        if len(_failed_stages) >= 3:
            warnings.warn(
                f"GSAS-II: {len(_failed_stages)} stages failed "
                f"({', '.join(_failed_stages)}). The refinement result may be "
                f"unreliable — check CIF quality or try simpler fitting range.")

        # ── Extract parameter uncertainties from covariance matrix ────────
        # After the final refinement, GSAS-II stores the covariance matrix
        # and per-parameter ESDs (estimated standard deviations) in the
        # project's Covariance data.  We extract scale-factor ESDs here
        # so we can propagate them into Hill & Howard weight fractions.
        #
        # GSAS-II parameter naming convention for phase scales:
        #   '<phase_id>::Scale'  (e.g. '0::Scale', '1::Scale')
        # The covariance data lives in gpx.data['Covariance']['data'].
        scale_sigmas = {}   # phase_name → sigma(Scale)
        try:
            # GSAS-II stores covariance data in multiple possible locations
            # depending on version.  Try them in order of preference.
            cov_data = {}
            for cov_path in [
                lambda: gpx.data.get('Covariance', {}).get('data', {}),
                lambda: gpx.data.get('Covariance', {}),
                lambda: gpx['Covariance']['data'] if 'Covariance' in gpx else {},
            ]:
                try:
                    _cd = cov_path()
                    if _cd and 'varyList' in _cd:
                        cov_data = _cd
                        break
                except (KeyError, TypeError, AttributeError):
                    continue

            if not cov_data:
                # Debug: log available top-level keys
                try:
                    top_keys = list(gpx.data.keys()) if hasattr(gpx.data, 'keys') else str(type(gpx.data))
                    print(f"  Covariance: not found. gpx.data keys: {top_keys}",
                          flush=True)
                except Exception:
                    print("  Covariance: not found, could not inspect gpx.data",
                          flush=True)

            vary_list = cov_data.get('varyList', [])
            sig_list  = cov_data.get('sig', [])
            # CRITICAL: sig_list is often a numpy array.  Do NOT use bare
            # 'if sig_list' — numpy arrays raise "ambiguous truth value"
            # when they contain >1 element.  Use len() instead.
            # Also convert to Python lists for safe .index() / 'in' usage.
            vary_list = list(vary_list)
            sig_list  = list(sig_list)
            if len(vary_list) > 0 and len(sig_list) > 0 and len(vary_list) == len(sig_list):
                # Log vary list for debugging scale parameter name format
                scale_vars = [v for v in vary_list if 'Scale' in v or 'scale' in v]
                print(f"  Covariance: {len(vary_list)} params, "
                      f"scale-related: {scale_vars}", flush=True)
                for idx_p, phase_obj in enumerate(gsas_phases):
                    # GSAS-II naming conventions for phase-histogram scales:
                    #   '{phase_id}:{hist_id}:Scale'  e.g. '0:0:Scale'
                    #   '{phase_id}::Scale'            e.g. '0::Scale'
                    # We try multiple patterns to be robust across versions.
                    found = False
                    for pattern in [f'{idx_p}:0:Scale', f'{idx_p}::Scale']:
                        if pattern in vary_list:
                            pos = vary_list.index(pattern)
                            scale_sigmas[phase_obj.name] = float(sig_list[pos])
                            found = True
                            break
                    if not found:
                        # Fallback: search for any Scale param for this phase
                        for vname, vsig in zip(vary_list, sig_list):
                            if vname.startswith(f'{idx_p}:') and 'Scale' in vname:
                                scale_sigmas[phase_obj.name] = float(vsig)
                                found = True
                                break
            else:
                print(f"  Covariance: varyList={len(vary_list) if vary_list else 0}, "
                      f"sig={len(sig_list) if sig_list else 0} — "
                      f"cannot extract ESDs", flush=True)
            if scale_sigmas:
                print(f"  Scale factor ESDs: {scale_sigmas}", flush=True)
            else:
                print("  Scale factor ESDs: not available (scale params "
                      "not found in vary list)", flush=True)
        except Exception as e:
            print(f"  Scale factor ESDs: extraction failed ({e})", flush=True)

        if progress_callback:
            progress_callback('GSAS-II: extracting results...')

        # ── Extract results ──────────────────────────────────────────────
        # Get calculated pattern
        y_calc_full = np.array(histogram.getdata('ycalc'))
        y_obs_full = np.array(histogram.getdata('yobs'))
        y_bg_full = np.array(histogram.getdata('background'))
        tt_full = np.array(histogram.getdata('x'))

        # Trim to our range
        rmask = (tt_full >= tt_min) & (tt_full <= tt_max)
        tt_out = tt_full[rmask]
        y_obs_out = y_obs_full[rmask]
        y_calc_out = y_calc_full[rmask]
        y_bg_out = y_bg_full[rmask]
        _y_bg_gsas = y_bg_out.copy()   # preserve for unbiased phase isolation

        # ── Background for display ──────────────────────────────────────
        # Use GSAS-II's actual refined background directly — no post-
        # processing.  A previous "dip correction" (trend + smoothing +
        # max) was applied here for display, but it created visual
        # artifacts (unphysical humps/discontinuities in the plotted
        # background) without affecting Rwp/GoF.  Disabled as of v5.
        y_bg_display = y_bg_out.copy()

        diff_out = y_obs_out - y_calc_out

        # ── REGRESSION CHECK: wtFactor must be included in chi²/GoF ────
        # Without this, reported GoF is inflated by K (sigma inflation
        # factor) even though Rwp/Rp are correct (overridden by GSAS-II's
        # own R-factors below).  This was a reporting bug found 2026-04-16.
        # If _SIGMA_INFLATION_FACTOR changes, this section automatically
        # stays in sync because we use the constant directly rather than
        # reading wtFactor back from histogram.data (GSAS-II may clear
        # or relocate that key during refinement cycles).
        base_weights = 1.0 / np.maximum(
            np.where(sig_r is not None, sig_r**2,
                     np.maximum(y_obs_out, 1.0)), 1e-6)
        weights_out = base_weights / (_SIGMA_INFLATION_FACTOR ** 2)

        # Compute statistics using GSAS-II's actual background
        n_params_est = sum(
            len(list(phase_obj.atoms())) + 7  # atoms + cell + scale + profile
            for phase_obj in gsas_phases
        ) + n_bg_coeffs + 1
        stats = compute_fit_statistics(y_obs_out, y_calc_out,
                                        weights_out, n_params_est)

        # Prefer GSAS-II's own R-factors (computed consistently with
        # its own background and weighting scheme)
        try:
            gsas_stats = histogram.get_statistics()
            if gsas_stats:
                stats['Rwp'] = round(gsas_stats.get('Rwp', stats['Rwp']), 2)
                stats['Rp'] = round(gsas_stats.get('Rp', stats['Rp']), 2)
        except Exception:
            pass

        # Use the display-corrected background for the output (plots)
        # but keep _y_bg_gsas for phase isolation (unchanged)
        y_bg_out = y_bg_display

        # Extract instrument params
        inst = _extract_instrument_params(histogram)

        # Zero shift
        try:
            zero_shift = float(
                histogram.data['Instrument Parameters'][0].get('Zero', [0, 0])[1])
        except Exception:
            zero_shift = 0.0

        # Sample displacement (refined Stage 2 + Stage 6).  For
        # bragg_brentano: DisplaceX (sample-height shift, mm).  For
        # capillary/Synergy-S: DisplaceY (radial offset from beam axis,
        # mm).  GSAS-II stores Sample Parameters as [value, refine_flag].
        # We expose this so users can confirm the cell-vs-displacement
        # trade-off didn't move position offset into the cell.
        try:
            _disp_param_name = ('DisplaceX' if geometry == 'bragg_brentano'
                                else 'DisplaceY')
            _disp_entry = histogram.data['Sample Parameters'].get(
                _disp_param_name, [0.0, False])
            displacement_um = float(_disp_entry[0]) if _disp_entry else 0.0
            displacement_param_name = _disp_param_name
        except Exception:
            displacement_um = 0.0
            displacement_param_name = 'unknown'

        # ── Per-phase results ────────────────────────────────────────────
        phase_patterns = []
        phase_results = []
        all_phase_refs = []   # per-phase reflection lists from generate_reflections

        # ── Weight fractions via Hill & Howard (1987) ────────────────────
        # W_α = S_α · Z_α · M_α · V_α  /  Σ(S_i · Z_i · M_i · V_i)
        # If Z or M are unavailable, fall back to raw scale normalisation.
        raw_scales = {}
        try:
            for phase_obj in gsas_phases:
                hapData = list(phase_obj.data['Histograms'].values())[0]
                scale_entry = hapData.get('Scale', [1.0])
                scale_val = scale_entry[0] if isinstance(scale_entry, (list, tuple)) else float(scale_entry)
                raw_scales[phase_obj.name] = scale_val
        except Exception as e:
            warnings.warn(f"GSAS-II: could not read scale factors: {e}")

        zmv_values = {}
        use_zmv = True
        for ph_input, phase_obj in zip(phases, gsas_phases):
            S = raw_scales.get(phase_obj.name, 1.0)
            _cell = phase_obj.data['General']['Cell']
            V_ph = cell_volume(float(_cell[1]), float(_cell[2]), float(_cell[3]),
                               float(_cell[4]), float(_cell[5]), float(_cell[6]))
            Z = ph_input.get('Z')
            M = molar_mass_from_formula(ph_input.get('formula', ''))
            if Z and M and V_ph > 0:
                zmv_values[phase_obj.name] = float(S) * float(Z) * float(M) * V_ph
            else:
                use_zmv = False
                break

        if not use_zmv:
            warnings.warn("GSAS-II: Z or molar mass unavailable for one or more "
                         "phases — falling back to raw scale factor weighting.")
            zmv_values = dict(raw_scales)

        total_zmv = sum(zmv_values.values()) or 1e-10
        print(f"GSAS-II scale factors: {raw_scales}", flush=True)
        print(f"GSAS-II ZMV values: {zmv_values}", flush=True)
        print(f"GSAS-II use_zmv: {use_zmv}", flush=True)

        # Warn if all scale factors are identical (common with failed refinement)
        scale_vals = list(raw_scales.values())
        if len(scale_vals) >= 2 and len(set(round(v, 6) for v in scale_vals)) == 1:
            warnings.warn("GSAS-II: all phase scale factors are identical — "
                         "refinement may not have converged properly.")

        # ── Weight fraction uncertainties via error propagation ──────────
        # For N phases with Hill & Howard weights:
        #   W_α = P_α / T   where P_α = S_α·Z_α·M_α·V_α, T = Σ P_i
        #
        # Partial derivative: ∂W_α/∂S_α = (Z_α·M_α·V_α) · (T - P_α) / T²
        #                     ∂W_α/∂S_β = -(Z_α·M_α·V_α·S_α) · (Z_β·M_β·V_β) / T²   (β≠α)
        #
        # FULL covariance treatment:
        # σ²(W_α) = Σ_i Σ_j (∂W_α/∂S_i) · (∂W_α/∂S_j) · Cov(S_i, S_j)
        #
        # For overlapping phases (WC + W2C), scale factors are strongly
        # ANTI-correlated (increasing one decreases the other).  Ignoring
        # this correlation drastically underestimates the uncertainties
        # (e.g. 0.05% instead of the true ~1-3%).
        wt_frac_sigmas = {}   # phase_name → σ(wt%) in percentage points
        wt_frac_sigma_sources = {}
        if scale_sigmas and len(zmv_values) >= 2:
            try:
                T = total_zmv
                phase_names = list(zmv_values.keys())
                n_ph = len(phase_names)

                # Build covariance matrix for scale factors from GSAS-II
                # If full covariance is available, use it; otherwise fall
                # back to diagonal (variances only) with a warning.
                cov_matrix = np.zeros((n_ph, n_ph))
                scale_indices = {}  # phase_name → index in vary_list
                have_full_cov = False
                try:
                    cov_data_mat = cov_data.get('covMatrix')
                    if cov_data_mat is not None:
                        cov_data_mat = np.array(cov_data_mat)
                        for idx_p, pname in enumerate(phase_names):
                            for pattern in [f'{idx_p}:0:Scale', f'{idx_p}::Scale']:
                                if pattern in vary_list:
                                    scale_indices[pname] = vary_list.index(pattern)
                                    break
                            if pname not in scale_indices:
                                for vname in vary_list:
                                    if vname.startswith(f'{idx_p}:') and 'Scale' in vname:
                                        scale_indices[pname] = vary_list.index(vname)
                                        break

                        if len(scale_indices) == n_ph:
                            for i_s, name_i in enumerate(phase_names):
                                for j_s, name_j in enumerate(phase_names):
                                    vi = scale_indices[name_i]
                                    vj = scale_indices[name_j]
                                    cov_matrix[i_s, j_s] = float(
                                        cov_data_mat[vi, vj])
                            have_full_cov = True
                            print(f"  Covariance: using full covariance matrix "
                                  f"for wt% uncertainties", flush=True)
                except Exception as e_cov:
                    print(f"  Covariance matrix extraction failed: {e_cov}",
                          flush=True)

                if not have_full_cov:
                    # Diagonal fallback — variances only
                    for i_s, pname in enumerate(phase_names):
                        sig_s = scale_sigmas.get(pname, 0.0)
                        cov_matrix[i_s, i_s] = sig_s ** 2
                    print(f"  Covariance: using diagonal only (no off-diagonal "
                          f"terms) — wt% uncertainties may be underestimated",
                          flush=True)

                # Compute ∂W_α/∂S_i for each alpha, then propagate
                for i_alpha, alpha in enumerate(phase_names):
                    P_alpha = zmv_values[alpha]
                    S_alpha = raw_scales.get(alpha, 1.0)
                    K_alpha = P_alpha / S_alpha if S_alpha > 0 else 0.0

                    # Build gradient vector dW_alpha/dS_i
                    grad = np.zeros(n_ph)
                    for i_beta, beta in enumerate(phase_names):
                        P_beta = zmv_values[beta]
                        S_beta = raw_scales.get(beta, 1.0)
                        K_beta = P_beta / S_beta if S_beta > 0 else 0.0
                        if beta == alpha:
                            grad[i_beta] = K_alpha * (T - P_alpha) / (T * T)
                        else:
                            grad[i_beta] = -P_alpha * K_beta / (T * T)

                    # σ²(W_α) = grad^T · Cov · grad
                    var_w = float(grad @ cov_matrix @ grad)
                    sigma_w_propagated = math.sqrt(max(var_w, 0.0)) * 100.0

                    # Systematic floor for routine lab XRD.
                    # Propagated ESD captures only the LM scale-factor
                    # uncertainty.  Real-sample wt% has additional
                    # systematic uncertainty from: (1) atomic Uiso held
                    # at default (no covariance), (2) sample mounting /
                    # absorption / preferred orientation residuals,
                    # (3) profile model imperfection (X/Y propagation
                    # to intensity isn't included), and (4) the choice
                    # between H&H and integrated-fraction reporting
                    # (the displayed wt% is integrated fraction, which
                    # is harder to bound rigorously).  Floor at the
                    # max of (propagated, 1.0% absolute, 2% relative).
                    _sys_floor = max(1.0, 0.02 * (P_alpha / T) * 100.0)
                    sigma_w = max(sigma_w_propagated, _sys_floor)
                    wt_frac_sigmas[alpha] = sigma_w
                    _cov_src = ('full_covariance' if have_full_cov
                                else 'diagonal_covariance')
                    if sigma_w > sigma_w_propagated * 1.01:
                        wt_frac_sigma_sources[alpha] = (
                            f'{_cov_src}+systematic_floor')
                    else:
                        wt_frac_sigma_sources[alpha] = _cov_src
                    if sigma_w > sigma_w_propagated * 1.01:
                        print(f"    {alpha}: propagated σ = "
                              f"±{sigma_w_propagated:.3f}%, "
                              f"reporting ±{sigma_w:.2f}% "
                              f"(systematic floor applied).",
                              flush=True)

                # For a closed binary mixture, wt%(B) = 100 - wt%(A), so
                # the absolute uncertainty in percentage points must be
                # the same for both displayed fractions. The per-phase
                # systematic floor above can otherwise make the major
                # phase report a larger floor than the minor phase.
                if n_ph == 2 and len(wt_frac_sigmas) == 2:
                    _binary_sigma = max(wt_frac_sigmas.values())
                    for _binary_name in phase_names:
                        wt_frac_sigmas[_binary_name] = _binary_sigma
                        _src = wt_frac_sigma_sources.get(
                            _binary_name, 'scale_covariance')
                        if 'binary_closure' not in _src:
                            wt_frac_sigma_sources[_binary_name] = (
                                f'{_src}+binary_closure')

                print(f"  Weight fraction uncertainties: "
                      f"{{ {', '.join(f'{k}: ±{v:.2f}%' for k, v in wt_frac_sigmas.items())} }}",
                      flush=True)
            except Exception as e:
                print(f"  Weight fraction uncertainty propagation failed: {e}",
                      flush=True)
        else:
            print("  Weight fraction uncertainties: not available "
                  "(need scale ESDs for ≥2 phases)", flush=True)

        # ── Extract GSAS-II RefList for physics-based profile generation ────
        # The RefList contains GSAS-II's refined Fc² values for each reflection,
        # giving more accurate phase envelopes than our generate_reflections output.
        # RefList columns: 0=h, 1=k, 2=l, 3=mult, 4=dsp, 5=2theta, 6=Fo²,
        #                  7=sig, 8=Fc², ...
        gsas_refs = {}   # phase_name → [(two_theta, d, (h,k,l), mult*Fc²)]
        try:
            raw_refl_lists = histogram.data.get('Reflection Lists', {})
            for ph_name, refl_data in raw_refl_lists.items():
                ref_arr = refl_data.get('RefList')
                if ref_arr is not None and len(ref_arr) > 0 and ref_arr.shape[1] > 8:
                    # First pass: collect all reflections and track max Fc²
                    # (raw, without multiplicity) for threshold filtering.
                    # Using raw Fc² matches how generate_reflections() filters
                    # in the preview — multiplicity should not determine whether
                    # a reflection is "real" vs "ghost".
                    raw_refs = []
                    max_fc2 = 0.0
                    for row in ref_arr:
                        h, k, l = int(row[0]), int(row[1]), int(row[2])
                        mult      = float(row[3])
                        d_sp      = float(row[4])
                        two_theta = float(row[5])
                        fc2       = float(row[8])   # Fc²
                        weight    = mult * fc2
                        if weight > 0 and tt_min <= two_theta <= tt_max:
                            raw_refs.append((two_theta, d_sp, (h, k, l), weight, fc2))
                            if fc2 > max_fc2:
                                max_fc2 = fc2
                    # Second pass: filter on raw Fc² (not mult*Fc²) to remove
                    # ghost reflections without losing legitimate low-multiplicity peaks.
                    # Use a very permissive threshold (0.01%) — heavy-atom carbide
                    # structures like W2C can have large Fc² dynamic range, and the
                    # previous 1% threshold was filtering out real reflections.
                    rel_thresh = max_fc2 * 1e-4  # 0.01% of strongest Fc²
                    refs = [(r[0], r[1], r[2], r[3])
                            for r in raw_refs if r[4] >= max(1e-6, rel_thresh)]
                    if refs:
                        gsas_refs[ph_name] = refs
                        print(f"  Loaded GSAS-II RefList for '{ph_name}': "
                              f"{len(refs)} reflections with Fc²", flush=True)
        except Exception as e:
            print(f"  (Could not extract GSAS-II RefList: {e})", flush=True)

        for i, (ph, phase_obj) in enumerate(zip(phases, gsas_phases)):
            # Cell parameters — read directly from GSAS-II data structure
            # to avoid get_cell() API differences across versions.
            # General['Cell'] = [mustRefine, a, b, c, alpha, beta, gamma, volume]
            _cell = phase_obj.data['General']['Cell']
            a, b, c = float(_cell[1]), float(_cell[2]), float(_cell[3])
            alpha, beta, gamma = float(_cell[4]), float(_cell[5]), float(_cell[6])

            # ── Sanity check: detect diverged cell parameters ────────────
            # If cell refinement failed and GSAS-II retained garbage values
            # (e.g. a=10⁸ Å from "Invalid cell metric tensor"), fall back
            # to the original input values from the phase dict.
            _orig_a = float(ph.get('a', a))
            _orig_b = float(ph.get('b', b))
            _orig_c = float(ph.get('c', c))
            _orig_alpha = float(ph.get('alpha', alpha))
            _orig_beta  = float(ph.get('beta', beta))
            _orig_gamma = float(ph.get('gamma', gamma))

            def _cell_diverged(ref, val, tol=2.0):
                """True if val deviates from ref by more than tol× factor."""
                if ref < 0.1:
                    return abs(val) > 100  # guard against near-zero ref
                ratio = val / ref
                return ratio > tol or ratio < (1.0 / tol)

            _any_diverged = (
                _cell_diverged(_orig_a, a) or
                _cell_diverged(_orig_b, b) or
                _cell_diverged(_orig_c, c) or
                abs(alpha - _orig_alpha) > 30 or
                abs(beta  - _orig_beta)  > 30 or
                abs(gamma - _orig_gamma) > 30
            )
            if _any_diverged:
                warnings.warn(
                    f"GSAS-II: cell params for phase {i} diverged "
                    f"(a={a:.4f} vs input {_orig_a:.4f}, "
                    f"b={b:.4f} vs {_orig_b:.4f}, "
                    f"c={c:.4f} vs {_orig_c:.4f}). "
                    f"Reverting to input values.")
                a, b, c = _orig_a, _orig_b, _orig_c
                alpha, beta, gamma = _orig_alpha, _orig_beta, _orig_gamma

            V = cell_volume(a, b, c, alpha, beta, gamma)
            V0 = cell_volume(_orig_a, _orig_b, _orig_c,
                             _orig_alpha, _orig_beta, _orig_gamma)

            def _pct_delta(val, ref):
                try:
                    val = float(val)
                    ref = float(ref)
                    if not math.isfinite(val) or not math.isfinite(ref):
                        return None
                    if abs(ref) < 1e-12:
                        return None
                    return (val - ref) / ref * 100.0
                except Exception:
                    return None

            _cell_delta = {
                'a_pct': _pct_delta(a, _orig_a),
                'b_pct': _pct_delta(b, _orig_b),
                'c_pct': _pct_delta(c, _orig_c),
                'volume_pct': _pct_delta(V, V0),
            }

            # Profile / size info
            prof = _extract_profile_params(phase_obj)
            cryst_A = prof.get('crystallite_size_A')
            cryst_source = 'gsas_hap_size' if cryst_A else None
            # In verification_mode, Stage 4b (size refinement) was skipped,
            # so HAP Size is at GSAS-II's default (~1 µm).  Discard that
            # value and let the Y-based fallback below estimate size from
            # the refined Y parameter.
            if verification_mode:
                cryst_A = None
                cryst_source = None

            # GSAS-II stores profile params in centidegrees² (sig) and
            # centidegrees (gam).  Convert to degrees² / degrees for
            # Caglioti/TCH compatibility with our functions.
            # 1 deg = 100 centideg  →  1 deg² = 10 000 centideg²
            U_deg = inst['U'] / 10000.0   # centideg² → deg²
            V_deg = inst['V'] / 10000.0
            W_deg = inst['W'] / 10000.0
            X_deg = inst['X'] / 100.0     # centideg → deg
            Y_deg = inst['Y'] / 100.0

            # If GSAS-II size extraction failed, estimate from Y
            if cryst_A is None and Y_deg > 0.01:
                cryst_A = size_from_Y(Y_deg, wavelength)
                if cryst_A:
                    cryst_source = 'Y_fallback'

            # FWHM and eta at a representative angle (40°)
            fwhm_reference_two_theta = 40.0
            fwhm_reference_hkl = None
            fwhm_reference_source = 'fallback_40deg'
            fwhm_rep, eta_rep = tch_fwhm_eta(
                fwhm_reference_two_theta, U_deg, V_deg, W_deg, X_deg, Y_deg)

            # Weight fraction (Hill & Howard or raw-scale fallback)
            scale_val = raw_scales.get(phase_obj.name, 1.0)
            zmv_val = zmv_values.get(phase_obj.name, scale_val)
            wt_pct = (zmv_val / total_zmv) * 100 if total_zmv > 0 else 0
            wt_pct_err = wt_frac_sigmas.get(phase_obj.name)  # may be None
            wt_pct_err_source = wt_frac_sigma_sources.get(phase_obj.name)

            # ── Generate tick positions / reflection list ─────────────────
            # Always use generate_reflections() for tick positions — this
            # matches the preview stick pattern and uses our pre-computed F²
            # with correct filtering.  GSAS-II's refined Fc² values can be
            # unreliable for phases that haven't fully converged (e.g. W2C
            # Pbcn), causing legitimate reflections to be filtered out.
            # The GSAS-II RefList (gsas_refs) is still used separately for
            # phase profile reconstruction in _compute_raw_phase_profile().
            #
            # CRITICAL: Use the GSAS-II CIF (gsas_cif_texts[i]) for site
            # expansion, NOT the original phase dict CIF.  When CifWriter
            # standardized the structure (e.g. W2C Pbcn), it may have
            # permuted axes.  The GSAS-II cell params (a,b,c above) are
            # in the CifWriter convention, so the sites must also be in
            # that convention for correct structure factor calculation.
            # Using the original CIF sites with CifWriter cell params →
            # wrong F(hkl) → real reflections filtered out → missing ticks.
            sys_ = (ph.get('system') or 'triclinic').lower()
            sg = ph.get('spacegroup_number', 1)
            _gsas_cif = gsas_cif_texts[i] if i < len(gsas_cif_texts) else ''
            sites = _get_expanded_sites(_gsas_cif or ph.get('cif_text', ''), sg)
            phase_refs = generate_reflections(
                a, b, c, alpha, beta, gamma, sys_, sg,
                wavelength, tt_min, tt_max, hkl_max=12,
                sites=sites)
            gsas_phase_refs = gsas_refs.get(phase_obj.name)

            # Report FWHM at a phase-relevant angle: the strongest
            # calculated reflection used for that phase's tick pattern.
            # The 40 deg value above is only a no-reflection fallback.
            try:
                _fwhm_ref = None
                if phase_refs:
                    _fwhm_ref = max(phase_refs, key=lambda r: float(r[3]))
                elif gsas_phase_refs:
                    _fwhm_ref = max(gsas_phase_refs,
                                    key=lambda r: float(r[3]))
                if _fwhm_ref:
                    fwhm_reference_two_theta = float(_fwhm_ref[0])
                    fwhm_reference_hkl = list(_fwhm_ref[2])
                    fwhm_reference_source = 'strongest_calculated_reflection'
                    fwhm_rep, eta_rep = tch_fwhm_eta(
                        fwhm_reference_two_theta,
                        U_deg, V_deg, W_deg, X_deg, Y_deg)
            except Exception:
                pass

            # Tick source: by default Python-filtered phase_refs (clean
            # display).  When use_gsas_ref_ticks is True, switch to
            # GSAS-II's RefList (all reflections including weak ones).
            if use_gsas_ticks_opt and gsas_phase_refs:
                tick_positions = [round(r[0], 3) for r in gsas_phase_refs
                                   if tt_min <= r[0] <= tt_max]
                _tick_src_label = 'GSAS-II RefList'
            else:
                tick_positions = [round(r[0], 3) for r in phase_refs]
                _tick_src_label = 'Python refs'
            print(f"  Tick positions for '{ph.get('name', '?')}' (SG {sg}): "
                  f"{len(tick_positions)} reflections in "
                  f"{tt_min:.1f}–{tt_max:.1f}° 2θ "
                  f"[source: {_tick_src_label}]", flush=True)
            if gsas_phase_refs:
                print(f"    GSAS-II RefList: {len(gsas_phase_refs)} "
                      f"reflections", flush=True)
            print(f"    Python refs:    {len(phase_refs)} reflections",
                  flush=True)
            # Manual reconstruction fallback uses Python refs.
            all_phase_refs.append(phase_refs)

            # B_iso (average over atoms; falls back to DEFAULT_B_ISO)
            b_iso_avg = DEFAULT_B_ISO
            try:
                atoms = list(phase_obj.atoms())
                if atoms:
                    uisos = [at.uiso for at in atoms if hasattr(at, 'uiso')]
                    if uisos:
                        b_iso_avg = 8 * math.pi**2 * np.mean(uisos)
            except Exception:
                pass

            # Preferred orientation is a reusable refinement output, not just
            # a submitted phase-card control.  Read the final GSAS-II HAP
            # value so the GUI can show the converged March-Dollase ratio for
            # the next run/preset.
            pref_ori_mode = 'off'
            pref_ori_value = None
            pref_ori_axis = None
            pref_ori_refined = False
            pref_ori_source = None
            try:
                hap_data = list(phase_obj.data['Histograms'].values())[0]
                pref_ori = hap_data.get('Pref.Ori.')
                if isinstance(pref_ori, (list, tuple)) and len(pref_ori) >= 2:
                    pref_ori_mode = _po_phase_modes.get(
                        i, 'global' if i in _pref_ori_phases else 'off')
                    pref_ori_value = float(pref_ori[1])
                    pref_ori_refined = (
                        bool(pref_ori[2]) if len(pref_ori) >= 3 else False)
                    if len(pref_ori) >= 4 and isinstance(
                            pref_ori[3], (list, tuple)):
                        pref_ori_axis = [int(x) for x in pref_ori[3]]
                    pref_ori_source = 'gsas_hap_pref_ori'
            except Exception:
                pass

            # Build a display name that always includes the space group symbol
            # so that phases with the same formula are clearly distinguishable.
            _sg_sym = (ph.get('spacegroup', '') or
                       _SG_HM.get(ph.get('spacegroup_number', 1), '') or
                       f"SG{ph.get('spacegroup_number', 1)}")
            _base_name = ph.get('name') or ph.get('formula') or f'Phase {i+1}'
            # Only append SG if not already present in the name
            if _sg_sym and _sg_sym.replace(' ', '') not in _base_name.replace(' ', ''):
                display_name = f"{_base_name} {_sg_sym}"
            else:
                display_name = _base_name

            phase_results.append({
                'name':              display_name,
                'cod_id':            ph.get('cod_id', ph.get('mp_id', '')),
                'formula':           ph.get('formula', ''),
                'a': round(a, 5), 'b': round(b, 5), 'c': round(c, 5),
                'alpha': round(alpha, 3), 'beta': round(beta, 3),
                'gamma': round(gamma, 3),
                'a0': round(_orig_a, 5), 'b0': round(_orig_b, 5),
                'c0': round(_orig_c, 5),
                'volume': round(V, 5),
                'volume0': round(V0, 5) if V0 else None,
                'delta_a_pct': (
                    round(_cell_delta['a_pct'], 3)
                    if _cell_delta['a_pct'] is not None else None),
                'delta_b_pct': (
                    round(_cell_delta['b_pct'], 3)
                    if _cell_delta['b_pct'] is not None else None),
                'delta_c_pct': (
                    round(_cell_delta['c_pct'], 3)
                    if _cell_delta['c_pct'] is not None else None),
                'delta_volume_pct': (
                    round(_cell_delta['volume_pct'], 3)
                    if _cell_delta['volume_pct'] is not None else None),
                'cell_reference': 'input_cif',
                'system':            (ph.get('system') or 'triclinic').lower(),
                'spacegroup_number': ph.get('spacegroup_number', 1),
                'spacegroup':        ph.get('spacegroup', ''),
                'scale':             round(scale_val, 5),
                'B_iso':             round(b_iso_avg, 4),
                'U': round(U_deg, 5), 'V': round(V_deg, 5),
                'W': round(W_deg, 5),
                'X': round(X_deg, 5), 'Y': round(Y_deg, 5),
                'preferred_orientation_mode': pref_ori_mode,
                'preferred_orientation_value': (
                    round(pref_ori_value, 5)
                    if pref_ori_value is not None else None),
                'preferred_orientation_axis': pref_ori_axis,
                'preferred_orientation_refined': pref_ori_refined,
                'preferred_orientation_source': pref_ori_source,
                'eta_at_strongest':  round(eta_rep, 3),
                'eta_at_40deg':      round(eta_rep, 3),
                'eta_at_fwhm_reference': round(eta_rep, 3),
                'fwhm_reference_two_theta': fwhm_reference_two_theta,
                'fwhm_reference_hkl': fwhm_reference_hkl,
                'fwhm_reference_source': fwhm_reference_source,
                'fwhm_deg':          round(fwhm_rep, 4),
                'crystallite_size_A':  round(cryst_A, 1) if cryst_A else None,
                'crystallite_size_nm': round(cryst_A / 10, 2) if cryst_A else None,
                'crystallite_size_source': cryst_source,
                'weight_fraction_%':       round(wt_pct, 1),
                'weight_fraction_err_%':   round(wt_pct_err, 2) if wt_pct_err is not None else None,
                'weight_fraction_err_source': wt_pct_err_source,
                'weight_fraction_method': 'hill_howard',
                'n_reflections':           len(tick_positions),
                'tick_positions':      tick_positions,
                'seeded_by':           'gsas2',
            })

        # ── Per-phase patterns ─────────────────────────────────────────
        # Primary method: GSAS-II phase isolation — zero all phases
        # except one, let GSAS-II recompute (0 cycles), extract
        # ycalc - background.  This uses GSAS-II's own profile
        # functions (correct peak shapes, asymmetry, everything).
        #
        # ── Phase decomposition ─────────────────────────────────────────
        # GSAS-II phase isolation is DISABLED for robustness.
        # It used GSAS-II's internal RefList for per-phase patterns,
        # which could disagree with the filtered Python reflection
        # set used for ticks.  Instead, we always use manual profile
        # reconstruction from all_phase_refs — same refs as ticks.
        # This guarantees phase envelopes match tick positions.

        if len(gsas_phases) > 1:
            total_above_bg = np.maximum(y_calc_out - y_bg_out, 0.0)
            decomp_ok = False

            # ── GSAS-II phase isolation: TOGGLEABLE ───────────────────
            # When ON: per-phase patterns come from GSAS-II's actual
            #   internal calculation (each phase's scale set to zero
            #   except one, ycalc recomputed, contribution extracted).
            #   These are the physically meaningful per-phase curves.
            # When OFF: per-phase patterns are reconstructed manually
            #   from the Python reflection list and refined U/V/W/X/Y.
            #   Display-only — guaranteed to match the tick positions
            #   but does not reflect GSAS-II's internal Fc² values.
            #
            # Rule: if isolation succeeds, USE IT.  If it fails (count
            # mismatch, exception, etc.), fall back to manual recon.
            # Never append both — the count guard at the end of the
            # block enforces this.
            _run_isolation = phase_isolation_opt
            if _run_isolation:
                print("  Phase isolation: ENABLED (per-phase curves from "
                      "GSAS-II ycalc).", flush=True)
            else:
                print("  Phase isolation: DISABLED (manual reconstruction "
                      "for tick/envelope consistency).", flush=True)

            # Save ALL refinement flags — the main refinement left many
            # params refinable (background, U/V/W/X/Y, cell, Uiso).
            # We must turn them ALL off before isolation, otherwise
            # do_refinements will try to refine them with one phase
            # zeroed, causing that phase's parameters to diverge.
            orig_hap_scales = []
            for phase_obj in gsas_phases:
                hapData = list(phase_obj.data['Histograms'].values())[0]
                orig_hap_scales.append(list(hapData['Scale']))

            # Save & turn off background refinement flag
            saved_bg_flag = histogram.data['Background'][0][1]
            histogram.data['Background'][0][1] = False

            # Save & turn off sample displacement refinement flags for
            # BOTH DisplaceX (Bragg-Brentano) and DisplaceY (Debye-Scherrer
            # transparency).  The pipeline only actively refines one of
            # these (DisplaceY for capillary data), but saving both is
            # defensive — if the geometry is changed later we don't want
            # a stale refinement flag to leak through phase isolation.
            saved_displace_flags = {}
            try:
                for _disp_key in ('DisplaceX', 'DisplaceY'):
                    _dp = histogram.data['Sample Parameters'].get(_disp_key)
                    if _dp and isinstance(_dp, list) and len(_dp) >= 2:
                        saved_displace_flags[_disp_key] = _dp[1]
                        _dp[1] = False
            except Exception:
                pass

            # Save & turn off instrument parameter refinement flags.
            # GSAS-II stores refine flags in TWO possible locations:
            # (a) inline: inst_params[0][key][2] (3-element lists)
            # (b) separate dict: inst_params[1][key] (boolean)
            inst_params_all = histogram.data['Instrument Parameters']
            inst_params_raw = inst_params_all[0]
            inst_refine_dict = (inst_params_all[1]
                                if len(inst_params_all) > 1
                                and isinstance(inst_params_all[1], dict)
                                else {})
            saved_inst_flags = {}
            for key in ['U', 'V', 'W', 'X', 'Y', 'SH/L', 'Zero']:
                # Location (a): inline in parameter list
                if key in inst_params_raw and len(inst_params_raw[key]) >= 3:
                    saved_inst_flags[(key, 'inline')] = \
                        inst_params_raw[key][2]
                    inst_params_raw[key][2] = False
                # Location (b): separate refine-flags dict
                if key in inst_refine_dict:
                    saved_inst_flags[(key, 'dict')] = \
                        inst_refine_dict[key]
                    inst_refine_dict[key] = False

            # Save & turn off preferred orientation refinement flags
            saved_pref_ori_flags = []
            for phase_obj in gsas_phases:
                try:
                    hapData = list(phase_obj.data['Histograms'].values())[0]
                    po = hapData.get('Pref.Ori.', ['MD', 1.0, False, [0,0,1]])
                    saved_pref_ori_flags.append(po[2])
                    po[2] = False
                except Exception:
                    saved_pref_ori_flags.append(False)


            # Save & turn off cell and atom refinement flags.
            # atom[2] = refinement flags string ('XU', 'U', '' etc.)
            # atom[9] = multiplicity — NOT the refinement flag.
            saved_cell_flags = []
            saved_atom_flags = []
            for phase_obj in gsas_phases:
                cell = phase_obj.data['General']['Cell']
                saved_cell_flags.append(cell[0])
                cell[0] = False
                atom_flags = []
                for atom in phase_obj.data['Atoms']:
                    if len(atom) > 2:
                        atom_flags.append(str(atom[2]))
                        atom[2] = ''   # clear XYZ + Uiso flags
                saved_atom_flags.append(atom_flags)

            # Belt-and-suspenders: force 0 cycles in Controls so GSAS-II
            # cannot iterate even if a flag was somehow missed.
            _controls = gpx.data.get('Controls', {}).get('data', {})
            _saved_maxcyc = _controls.get('max cyc', 3)
            _controls['max cyc'] = 0

            try:
                phase_patterns = []
                if not _run_isolation:
                    raise _IsolationSkipped()
                for i in range(len(gsas_phases)):
                    # Activate only phase i; zero all others
                    for j, phase_obj in enumerate(gsas_phases):
                        hapData = list(phase_obj.data['Histograms'].values())[0]
                        if j == i:
                            hapData['Scale'] = [orig_hap_scales[j][0], False]
                        else:
                            hapData['Scale'] = [0.0, False]

                    # Persist to disk so do_refinements picks up changes.
                    gpx.save()

                    # Evaluate only (cycles=0) — all refinement flags are
                    # OFF, so nothing can diverge.
                    gpx.do_refinements([{'set': {}, 'cycles': 0}])

                    # Re-fetch histogram — do_refinements reloads the
                    # project from disk, so the old reference may be stale.
                    histogram = gpx.histograms()[0]

                    y_calc_iso = np.array(histogram.getdata('ycalc'))
                    # Use GSAS-II's original background (before our dip
                    # correction) so the correction delta doesn't bias
                    # phase ratios toward the dominant phase.
                    phase_pat = np.maximum(
                        y_calc_iso[rmask] - _y_bg_gsas, 0.0)
                    phase_patterns.append(phase_pat.tolist())

                    ph_name = (phase_results[i]['name']
                               if i < len(phase_results) else f'Phase {i}')
                    print(f"  Phase {i} ({ph_name}): "
                          f"max={np.max(phase_pat):.1f}, "
                          f"integrated={np.sum(phase_pat):.1f}", flush=True)

                # Diagnostic: verify isolation produced correct patterns.
                sum_iso = np.zeros_like(tt_out, dtype=np.float64)
                for pp in phase_patterns:
                    sum_iso += np.array(pp)
                max_diff = float(np.max(np.abs(
                    sum_iso - total_above_bg)))
                sum_total = float(np.sum(total_above_bg))
                sum_iso_total = float(np.sum(sum_iso))
                print(f"  Diagnostic: sum(iso)={sum_iso_total:.1f}, "
                      f"sum(above_bg)={sum_total:.1f}, "
                      f"max_pointwise_diff={max_diff:.2f}", flush=True)

                # Proportional normalization: ensure stacked fills sum
                # exactly to total_above_bg at every point so that the
                # fills match the I_calc line perfectly.
                for i_pp, pp in enumerate(phase_patterns):
                    pp_arr = np.array(pp)
                    ratio = np.zeros_like(sum_iso)
                    np.divide(pp_arr, sum_iso, out=ratio,
                              where=sum_iso > 0)
                    phase_patterns[i_pp] = np.maximum(
                        ratio * total_above_bg, 0.0).tolist()

                decomp_ok = True
                print("  Phase isolation succeeded.", flush=True)
            except _IsolationSkipped:
                phase_patterns = []  # silently skip to fallback
            except Exception as e_iso:
                print(f"  Phase isolation failed: {e_iso}", flush=True)
                phase_patterns = []
            finally:
                # Always restore ALL original flags and scales
                histogram.data['Background'][0][1] = saved_bg_flag
                for _disp_key, _disp_flag in saved_displace_flags.items():
                    try:
                        histogram.data['Sample Parameters'][_disp_key][1] = \
                            _disp_flag
                    except Exception:
                        pass
                for (key, loc), flag in saved_inst_flags.items():
                    if loc == 'inline':
                        inst_params_raw[key][2] = flag
                    else:
                        inst_refine_dict[key] = flag
                for idx_r, phase_obj in enumerate(gsas_phases):
                    if idx_r < len(saved_cell_flags):
                        phase_obj.data['General']['Cell'][0] = \
                            saved_cell_flags[idx_r]
                    for j_a, atom in enumerate(phase_obj.data['Atoms']):
                        if (idx_r < len(saved_atom_flags)
                                and len(atom) > 2
                                and j_a < len(saved_atom_flags[idx_r])):
                            atom[2] = saved_atom_flags[idx_r][j_a]
                    hapData = list(phase_obj.data['Histograms'].values())[0]
                    if idx_r < len(orig_hap_scales):
                        hapData['Scale'] = orig_hap_scales[idx_r]
                    # Restore preferred orientation refinement flag
                    if idx_r < len(saved_pref_ori_flags):
                        try:
                            po = hapData.get('Pref.Ori.')
                            if po:
                                po[2] = saved_pref_ori_flags[idx_r]
                        except Exception:
                            pass
                # Restore max cycles to original value
                _controls['max cyc'] = _saved_maxcyc
                gpx.save()
                histogram = gpx.histograms()[0]

            # ── Fallback: profile reconstruction from reflections ────
            if not decomp_ok:
                print("  Falling back to profile reconstruction...",
                      flush=True)
                if all_phase_refs and len(all_phase_refs) == len(gsas_phases):
                    try:
                        U_d = inst['U'] / 10000.0
                        V_d = inst['V'] / 10000.0
                        W_d = inst['W'] / 10000.0
                        X_d = inst['X'] / 100.0
                        Y_d = inst['Y'] / 100.0

                        raw_profiles = []
                        for i_ph, (phase_obj_r, fallback_refs) in enumerate(
                                zip(gsas_phases, all_phase_refs)):
                            # Always use all_phase_refs (same as ticks)
                            # for consistency.  Never mix with GSAS-II
                            # RefList which may have different reflections.
                            refs_to_use = fallback_refs
                            if not refs_to_use:
                                raise ValueError(
                                    f"No reflections for phase {i_ph}")
                            raw_profiles.append(
                                _compute_raw_phase_profile(
                                    tt_out, refs_to_use,
                                    U_d, V_d, W_d, X_d, Y_d,
                                    gaussian_only=True))

                        scaled = []
                        for phase_obj, prof in zip(gsas_phases, raw_profiles):
                            s = raw_scales.get(phase_obj.name, 1.0)
                            scaled.append(prof * s)

                        sum_raw = np.zeros_like(tt_out, dtype=np.float64)
                        for sp in scaled:
                            sum_raw += sp
                        peak_max = (np.max(sum_raw)
                                    if np.max(sum_raw) > 0 else 1.0)
                        threshold = peak_max * 1e-10
                        for sp in scaled:
                            ratio = np.zeros_like(sum_raw)
                            np.divide(sp, sum_raw, out=ratio,
                                      where=sum_raw > threshold)
                            phase_patterns.append(
                                np.maximum(
                                    ratio * total_above_bg, 0.0).tolist())
                        decomp_ok = True
                        print("  Profile reconstruction succeeded.",
                              flush=True)
                    except Exception as e_fc2:
                        print(f"  Profile reconstruction failed: {e_fc2}",
                              flush=True)
                        phase_patterns = []

            # ── Last resort: equal split ───────────────────────────────
            if not decomp_ok:
                warnings.warn(
                    "GSAS-II phase decomposition: both isolation and "
                    "profile reconstruction failed. Falling back to "
                    "equal split — weight fractions will be UNRELIABLE.")
                phase_patterns = []
                for _ in gsas_phases:
                    share = (total_above_bg / len(gsas_phases)).tolist()
                    phase_patterns.append(share)

            # Hard guard: phase_patterns must match gsas_phases count
            if len(phase_patterns) != len(gsas_phases):
                print(f"  WARNING: phase pattern count mismatch: "
                      f"{len(phase_patterns)} patterns for "
                      f"{len(gsas_phases)} phases. Discarding "
                      f"decomposition.", flush=True)
                phase_patterns = []
                decomp_ok = False

            # Use integrated pattern fractions as the displayed weight
            # fractions.  The Hill & Howard (H&H) method is the standard
            # Rietveld approach, but it requires accurate Z (formula units
            # per cell) and molar mass for each phase.  When these are
            # incorrect or when scale factors haven't fully converged,
            # H&H can give wildly wrong values.  The integrated fractions
            # are a direct measure of each phase's contribution to the
            # diffraction signal and are more robust for display.
            #
            # NOTE: the wt% UNCERTAINTIES are still derived from H&H
            # error propagation (covariance matrix).  These are logged
            # alongside the H&H values for diagnostics.
            total_integ = sum(
                np.sum(np.array(pp)) for pp in phase_patterns) or 1e-30
            for i_wp, pp in enumerate(phase_patterns):
                integ_frac = np.sum(np.array(pp)) / total_integ * 100
                hh_frac = (phase_results[i_wp]['weight_fraction_%']
                           if i_wp < len(phase_results) else None)
                print(f"  Phase {i_wp}: "
                      f"H&H wt% = {hh_frac}%, "
                      f"integrated fraction = {integ_frac:.1f}%", flush=True)
                # Overwrite with integrated fraction for display
                if i_wp < len(phase_results):
                    phase_results[i_wp]['weight_fraction_%'] = round(integ_frac, 1)
                    phase_results[i_wp]['weight_fraction_method'] = (
                        'integrated_phase_pattern')
        else:
            # Single phase — entire signal above background
            total_above_bg = np.maximum(y_calc_out - y_bg_out, 0.0)
            phase_patterns.append(total_above_bg.tolist())

        # Use the actual per-phase calculated pattern to choose the FWHM
        # reference angle shown in the GUI. This avoids reporting a generic
        # 40 deg width when the phase pattern has a clear strongest peak.
        try:
            _tt_for_fwhm = np.asarray(tt_out, dtype=float)
            for _pi, _pp in enumerate(phase_patterns):
                if _pi >= len(phase_results):
                    continue
                _pp_arr = np.asarray(_pp, dtype=float)
                _n = min(_pp_arr.size, _tt_for_fwhm.size)
                if _n <= 0:
                    continue
                _pk_idx = int(np.argmax(_pp_arr[:_n]))
                if (not np.isfinite(_pp_arr[_pk_idx])
                        or float(_pp_arr[_pk_idx]) <= 0):
                    continue
                _tt_ref = float(_tt_for_fwhm[_pk_idx])
                _pr = phase_results[_pi]
                _U = float(_pr.get('U', 0.0) or 0.0)
                _V = float(_pr.get('V', 0.0) or 0.0)
                _W = float(_pr.get('W', 0.0) or 0.0)
                _X = float(_pr.get('X', 0.0) or 0.0)
                _Y = float(_pr.get('Y', 0.0) or 0.0)
                _fwhm, _eta = tch_fwhm_eta(_tt_ref, _U, _V, _W, _X, _Y)
                _pr['fwhm_reference_two_theta'] = round(_tt_ref, 4)
                _pr['fwhm_reference_source'] = (
                    'strongest_phase_pattern_peak')
                _pr['fwhm_deg'] = round(_fwhm, 4)
                _pr['eta_at_fwhm_reference'] = round(_eta, 3)
                _refs = all_phase_refs[_pi] if _pi < len(all_phase_refs) else []
                if _refs:
                    _nearest = min(_refs, key=lambda r: abs(float(r[0]) - _tt_ref))
                    if abs(float(_nearest[0]) - _tt_ref) <= 0.5:
                        _pr['fwhm_reference_hkl'] = list(_nearest[2])
        except Exception as _e:
            print(f"  FWHM reference update from phase patterns failed: {_e}",
                  flush=True)

        # ── Per-phase Scherrer fallback ────────────────────────────────
        # When HAP Size wasn't refined for a phase, estimate crystallite
        # size from an isolated peak in that phase's contribution using
        # the Scherrer equation:  D = K·λ / (β·cos θ),  K = 0.9.
        #
        # Picks the strongest peak in the phase pattern and validates
        # that other phases' contributions at that 2θ are small (< 25%
        # of this phase's value), so the FWHM measurement isn't
        # contaminated by overlap.
        try:
            _SCHERRER_K = 0.9
            _PI = math.pi
            tt_arr = np.asarray(tt_out, dtype=float)
            _phase_patterns_arr = [np.asarray(pp, dtype=float)
                                    for pp in phase_patterns]
            for _pi, _pp in enumerate(_phase_patterns_arr):
                if _pi >= len(phase_results):
                    continue
                # Skip if size is already populated from refinement / Y-fallback
                if phase_results[_pi].get('crystallite_size_A'):
                    continue
                if _pp.size < 3 or float(np.max(_pp)) <= 0:
                    continue
                # Find candidate peak: strongest local max with FWHM ≥ 2 pts
                # AND other phases' contribution there < 25% of this phase
                _pk_max = float(np.max(_pp))
                _pk_idx = int(np.argmax(_pp))
                # Overlap check: sum of OTHER phase patterns at peak position
                _other_sum = 0.0
                for _oi, _opp in enumerate(_phase_patterns_arr):
                    if _oi == _pi:
                        continue
                    if _pk_idx < len(_opp):
                        _other_sum += float(_opp[_pk_idx])
                _overlap_frac = _other_sum / max(_pk_max, 1e-12)
                if _overlap_frac > 0.25:
                    print(f"  Scherrer (phase {_pi}): strongest peak at "
                          f"2θ={float(tt_arr[_pk_idx]):.2f}° has "
                          f"{_overlap_frac*100:.0f}% overlap with other "
                          f"phases — skipping.", flush=True)
                    continue
                # Measure FWHM around _pk_idx
                _half = _pk_max / 2.0
                # Walk left
                _li = _pk_idx
                while _li > 0 and _pp[_li] > _half:
                    _li -= 1
                # Walk right
                _ri = _pk_idx
                while _ri < len(_pp) - 1 and _pp[_ri] > _half:
                    _ri += 1
                if _ri <= _li + 1:
                    continue
                _fwhm_deg = float(tt_arr[_ri] - tt_arr[_li])
                _tt_peak  = float(tt_arr[_pk_idx])
                if _fwhm_deg <= 0 or _tt_peak <= 0:
                    continue
                # Scherrer.  Subtract instrumental broadening if we know
                # it from the measured instprm Y (centideg) — gives a
                # sample-only FWHM.  Conservative: use the full FWHM if
                # instrumental contribution is unknown.
                _theta_rad = math.radians(_tt_peak / 2.0)
                _beta_rad  = math.radians(_fwhm_deg)
                _D_A = _SCHERRER_K * wavelength / (_beta_rad
                                                    * math.cos(_theta_rad))
                if _D_A > 0 and _D_A < 1e6:
                    phase_results[_pi]['crystallite_size_A'] = round(
                        float(_D_A), 1)
                    phase_results[_pi]['crystallite_size_nm'] = round(
                        float(_D_A) / 10.0, 2)
                    phase_results[_pi]['size_method'] = 'Scherrer'
                    phase_results[_pi]['crystallite_size_source'] = (
                        'scherrer_phase_pattern')
                    print(f"  Scherrer (phase {_pi}): peak at "
                          f"2θ={_tt_peak:.2f}°, FWHM={_fwhm_deg:.3f}°, "
                          f"D = {_D_A:.0f} Å = {_D_A/10:.1f} nm "
                          f"(overlap {_overlap_frac*100:.1f}%).",
                          flush=True)
        except Exception as _e:
            print(f"  Scherrer fallback failed: {_e}", flush=True)

        # ── Post-refinement sanity warnings ─────────────────────────────
        # Flag suspicious results without failing the run.  Useful for
        # automated screening where a low Rwp is not always trustworthy.
        _sanity_warnings = []

        # Zero shift
        if _y_nonnegative_clamped:
            _sanity_warnings.append(
                "Refined Y became negative and was clamped to 0; this is a "
                "sign that X/Y broadening terms may be trading off.")

        # High-leverage knobs enabled in a new-sample recipe
        if refine_xyz:
            _sanity_warnings.append(
                "XYZ atom-position refinement was enabled. Lab XRD usually "
                "does not constrain atom coordinates strongly enough unless "
                "the structure, instrument profile, and constraints are "
                "already well validated; XYZ can otherwise absorb CIF or "
                "instrument errors and make a wrong structure look better.")
        if refine_uiso_opt:
            _sanity_warnings.append(
                "Uiso refinement was enabled. Uiso changes calculated peak "
                "intensities through the Debye-Waller factor, especially at "
                "higher angle; a lower Rwp can therefore come from intensity "
                "reweighting rather than a better phase fraction. Check that "
                "B_iso values are physically plausible and that wt% changes "
                "remain stable against the baseline.")
        if refine_size_opt:
            _sanity_warnings.append(
                "Global size refinement was enabled. Size broadening can "
                "compensate for instrument/profile errors or unresolved "
                "microstrain, so keep it only if peak widths improve in a "
                "physically expected way versus the baseline.")
        for _idx_opt, _popts_warn in enumerate(phase_options_list or []):
            if not isinstance(_popts_warn, dict):
                continue
            if _popts_warn.get('refine_size'):
                _sanity_warnings.append(
                    f"Phase {_idx_opt + 1} size refinement is active. It "
                    f"changes that phase's peak widths and can trade off "
                    f"with profile/instrument broadening; confirm this phase "
                    f"is visibly broader than the others.")
            if _popts_warn.get('refine_mustrain'):
                _sanity_warnings.append(
                    f"Phase {_idx_opt + 1} microstrain refinement is active. "
                    f"Microstrain and size both broaden peaks, but with "
                    f"different angle/hkl trends; confirm the residuals show "
                    f"strain-like broadening before interpreting it.")
            if _popts_warn.get('po_mode') == 'refined':
                _sanity_warnings.append(
                    f"Phase {_idx_opt + 1} preferred orientation is refined. "
                    f"PO changes relative reflection intensities and can "
                    f"trade off with phase fraction; compare against off or "
                    f"fixed-PO fits and check wt% stability.")
            if _popts_warn.get('uniform_cell'):
                _sanity_warnings.append(
                    f"Phase {_idx_opt + 1} uniform-cell post-scaling is "
                    f"active. This is a diagnostic constraint: small fit "
                    f"penalty supports uniform contraction/expansion, while "
                    f"a large penalty suggests real anisotropic cell change "
                    f"or an inadequate model.")

        for pr in phase_results:
            _name_warn = pr.get('name', '?')
            for _axis_label, _key_warn in (
                    ('a', 'delta_a_pct'), ('b', 'delta_b_pct'),
                    ('c', 'delta_c_pct'), ('V', 'delta_volume_pct')):
                _dv = pr.get(_key_warn)
                if _dv is not None and abs(float(_dv)) > 3.0:
                    _sanity_warnings.append(
                        f"Large {_axis_label} cell change for '{_name_warn}' "
                        f"({_dv}% vs input CIF); check phase identity and "
                        f"sample displacement/zero correction.")
            _wf = pr.get('weight_fraction_%')
            _we = pr.get('weight_fraction_err_%')
            if _wf is not None and _we is not None:
                if float(_we) > 5.0 or (
                        abs(float(_wf)) > 1e-9
                        and float(_we) / abs(float(_wf)) > 0.25):
                    _sanity_warnings.append(
                        f"Large weight-fraction uncertainty for "
                        f"'{_name_warn}' ({_wf} +/- {_we} wt%).")
            _src_warn = str(pr.get('crystallite_size_source') or '').lower()
            if 'fallback' in _src_warn:
                _sanity_warnings.append(
                    f"Crystallite size for '{_name_warn}' came from a "
                    f"fallback estimate; do not treat it as a direct "
                    f"refined size parameter.")

        # Zero shift
        if abs(zero_shift) > 0.1:
            _sanity_warnings.append(
                f"Large zero shift ({zero_shift:.4f}°); check wavelength "
                f"or instrument geometry.")

        # Displacement
        try:
            _disp_val = float(histogram.data['Sample Parameters']
                              .get(_displace_param, [0])[0])
            if abs(_disp_val) > 500:
                _sanity_warnings.append(
                    f"Large displacement correction ({_displace_param}="
                    f"{_disp_val:.1f} µm); check {geometry} geometry.")
        except Exception:
            pass

        # March-Dollase ratio
        for idx, phase_obj in enumerate(gsas_phases):
            try:
                hapData = list(phase_obj.data['Histograms'].values())[0]
                md_ratio = hapData.get('Pref.Ori.', ['MD', 1.0])[1]
                if md_ratio < 0.5 or md_ratio > 2.0:
                    _sanity_warnings.append(
                        f"Extreme March-Dollase ratio ({md_ratio:.3f}) for "
                        f"phase {idx}; weight fractions may be unreliable.")
            except Exception:
                pass

        # Near-zero weight fraction
        for pr in phase_results:
            wf = pr.get('weight_fraction_%', 50)
            if wf is not None and wf < 1.0:
                _sanity_warnings.append(
                    f"Phase '{pr.get('name', '?')}' refined to {wf:.1f} wt% "
                    f"— consider removing it from the model.")

        if _sanity_warnings:
            for _sw in _sanity_warnings:
                print(f"  WARNING: {_sw}", flush=True)

        # ── Hard validation: one source of truth ─────────────────────
        _n_phases = len(gsas_phases)
        if len(phase_patterns) != _n_phases:
            print(f"  WARNING: phase_patterns ({len(phase_patterns)}) "
                  f"!= gsas_phases ({_n_phases})", flush=True)
        if len(phase_results) != _n_phases:
            print(f"  WARNING: phase_results ({len(phase_results)}) "
                  f"!= gsas_phases ({_n_phases})", flush=True)

        result = {
            'tt':             tt_out.tolist(),
            'y_obs':          y_obs_out.tolist(),
            'y_calc':         y_calc_out.tolist(),
            'y_background':   y_bg_out.tolist(),
            'phase_patterns': phase_patterns,
            'residuals':      diff_out.tolist(),
            'statistics':     stats,
            'phase_results':  phase_results,
            'zero_shift':     round(zero_shift, 5),
            'displacement_um': round(displacement_um, 2),
            'displacement_param': displacement_param_name,
            'wavelength':     wavelength,
            'pymatgen_used':  False,
            'method':         'GSAS-II',
            'warnings':       _sanity_warnings,
            'instrument':     instrument,
            'instrument_label': profile['label'],
            'instrument_reason': instrument_reason,
        }

        print("=== GSAS-II REFINEMENT DONE ===", flush=True)
        return result

    finally:
        # Clean up temporary files
        for p in cif_paths:
            try:
                os.unlink(p)
            except OSError:
                pass
        try:
            os.unlink(data_path)
        except OSError:
            pass
        # Clean up temp working directory (including .gpx, .bak files)
        try:
            import shutil
            shutil.rmtree(work_dir, ignore_errors=True)
        except Exception:
            pass
