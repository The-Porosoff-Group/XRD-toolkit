#!/usr/bin/env python
"""
calibrate_instprm.py — Generate a measured .instprm from NIST SRM 640g Si

Uses GSAS-II's own profile model to fit the Si standard, then writes a
.instprm file with the refined U/V/W/X/Y/SH/L/Zero parameters.

Run from C:\catalysis-toolkit:
    .conda_env\python.exe calibrate_instprm.py
"""

import os, sys, tempfile, shutil, traceback
import numpy as np

# ── Paths (relative to THIS script's location) ────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

SI_XY_PATH = os.path.join(SCRIPT_DIR,
    'POROSOFF_NIST-Si640',
    'BB_0.01st-12sp-0.5IS-20RS1-OpenRS2--NIST-Si640_001.xy')

OUTPUT_INSTPRM = os.path.join(SCRIPT_DIR, 'smartlab_Si640g.instprm')

# NIST SRM 640g Si: Fd-3m (#227), a = 5.431109 A
SI_A = 5.431109
LAM1 = 1.540593
LAM2 = 1.544414
I_RATIO = 0.5
POLARIZ = 0.5


def log(msg):
    print(msg, flush=True)


def main():
    log(f"Script dir: {SCRIPT_DIR}")

    # ── Find data file ─────────────────────────────────────────────────
    if not os.path.isfile(SI_XY_PATH):
        log(f"ERROR: Si data not found at:\n  {SI_XY_PATH}")
        log("Copy the POROSOFF_NIST-Si640 folder into C:\\catalysis-toolkit\\")
        sys.exit(1)
    log(f"Si 640g data: {SI_XY_PATH}")
    log(f"Output:       {OUTPUT_INSTPRM}")

    # ── Load data ──────────────────────────────────────────────────────
    data = np.loadtxt(SI_XY_PATH)
    tt = data[:, 0]
    intensity = data[:, 1]
    sigma = np.sqrt(np.maximum(intensity, 1.0))
    log(f"Data: {len(tt)} points, {tt[0]:.2f} - {tt[-1]:.2f} deg")

    # ── Set up temp directory ──────────────────────────────────────────
    work_dir = tempfile.mkdtemp(prefix='si640g_cal_')
    log(f"Working dir: {work_dir}")

    # Write .xye
    xye_path = os.path.join(work_dir, 'si640g.xye')
    with open(xye_path, 'w') as f:
        for i in range(len(tt)):
            f.write(f"{tt[i]:.6f}  {intensity[i]:.4f}  {sigma[i]:.4f}\n")
    log("Wrote .xye file")

    # Write Si CIF
    cif_path = os.path.join(work_dir, 'Si_640g.cif')
    with open(cif_path, 'w') as f:
        f.write(f"""data_Si_640g
_chemical_formula_sum 'Si'
_cell_formula_units_Z 8
_cell_length_a {SI_A}
_cell_length_b {SI_A}
_cell_length_c {SI_A}
_cell_angle_alpha 90.000
_cell_angle_beta 90.000
_cell_angle_gamma 90.000
_symmetry_Int_Tables_number 227
_symmetry_space_group_name_H-M 'F d -3 m'
loop_
_symmetry_equiv_pos_as_xyz
'x, y, z'
'-x, -y, z'
'-x, y, -z'
'x, -y, -z'
'z, x, y'
'z, -x, -y'
'-z, -x, y'
'-z, x, -y'
'y, z, x'
'-y, z, -x'
'y, -z, -x'
'-y, -z, x'
'-y, -x, -z'
'y, x, -z'
'-y, x, z'
'y, -x, z'
'-x, -z, -y'
'x, z, -y'
'x, -z, y'
'-x, z, y'
'-z, -y, -x'
'-z, y, x'
'z, y, -x'
'z, -y, x'
'x+1/2, y+1/2, z'
'-x+1/2, -y+1/2, z'
'-x+1/2, y+1/2, -z'
'x+1/2, -y+1/2, -z'
'z+1/2, x+1/2, y'
'z+1/2, -x+1/2, -y'
'-z+1/2, -x+1/2, y'
'-z+1/2, x+1/2, -y'
'y+1/2, z+1/2, x'
'-y+1/2, z+1/2, -x'
'y+1/2, -z+1/2, -x'
'-y+1/2, -z+1/2, x'
'-y+1/2, -x+1/2, -z'
'y+1/2, x+1/2, -z'
'-y+1/2, x+1/2, z'
'y+1/2, -x+1/2, z'
'-x+1/2, -z+1/2, -y'
'x+1/2, z+1/2, -y'
'x+1/2, -z+1/2, y'
'-x+1/2, z+1/2, y'
'-z+1/2, -y+1/2, -x'
'-z+1/2, y+1/2, x'
'z+1/2, y+1/2, -x'
'z+1/2, -y+1/2, x'
'x+1/2, y, z+1/2'
'-x+1/2, -y, z+1/2'
'-x+1/2, y, -z+1/2'
'x+1/2, -y, -z+1/2'
'z+1/2, x, y+1/2'
'z+1/2, -x, -y+1/2'
'-z+1/2, -x, y+1/2'
'-z+1/2, x, -y+1/2'
'y+1/2, z, x+1/2'
'-y+1/2, z, -x+1/2'
'y+1/2, -z, -x+1/2'
'-y+1/2, -z, x+1/2'
'-y+1/2, -x, -z+1/2'
'y+1/2, x, -z+1/2'
'-y+1/2, x, z+1/2'
'y+1/2, -x, z+1/2'
'-x+1/2, -z, -y+1/2'
'x+1/2, z, -y+1/2'
'x+1/2, -z, y+1/2'
'-x+1/2, z, y+1/2'
'-z+1/2, -y, -x+1/2'
'-z+1/2, y, x+1/2'
'z+1/2, y, -x+1/2'
'z+1/2, -y, x+1/2'
'x, y+1/2, z+1/2'
'-x, -y+1/2, z+1/2'
'-x, y+1/2, -z+1/2'
'x, -y+1/2, -z+1/2'
'z, x+1/2, y+1/2'
'z, -x+1/2, -y+1/2'
'-z, -x+1/2, y+1/2'
'-z, x+1/2, -y+1/2'
'y, z+1/2, x+1/2'
'-y, z+1/2, -x+1/2'
'y, -z+1/2, -x+1/2'
'-y, -z+1/2, x+1/2'
'-y, -x+1/2, -z+1/2'
'y, x+1/2, -z+1/2'
'-y, x+1/2, z+1/2'
'y, -x+1/2, z+1/2'
'-x, -z+1/2, -y+1/2'
'x, z+1/2, -y+1/2'
'x, -z+1/2, y+1/2'
'-x, z+1/2, y+1/2'
'-z, -y+1/2, -x+1/2'
'-z, y+1/2, x+1/2'
'z, y+1/2, -x+1/2'
'z, -y+1/2, x+1/2'
'x+1/2, y+1/2, z+1/2'
'-x+1/2, -y+1/2, z+1/2'
'-x+1/2, y+1/2, -z+1/2'
'x+1/2, -y+1/2, -z+1/2'
'z+1/2, x+1/2, y+1/2'
'z+1/2, -x+1/2, -y+1/2'
'-z+1/2, -x+1/2, y+1/2'
'-z+1/2, x+1/2, -y+1/2'
'y+1/2, z+1/2, x+1/2'
'-y+1/2, z+1/2, -x+1/2'
'y+1/2, -z+1/2, -x+1/2'
'-y+1/2, -z+1/2, x+1/2'
'-y+1/2, -x+1/2, -z+1/2'
'y+1/2, x+1/2, -z+1/2'
'-y+1/2, x+1/2, z+1/2'
'y+1/2, -x+1/2, z+1/2'
'-x+1/2, -z+1/2, -y+1/2'
'x+1/2, z+1/2, -y+1/2'
'x+1/2, -z+1/2, y+1/2'
'-x+1/2, z+1/2, y+1/2'
'-z+1/2, -y+1/2, -x+1/2'
'-z+1/2, y+1/2, x+1/2'
'z+1/2, y+1/2, -x+1/2'
'z+1/2, -y+1/2, x+1/2'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
Si1  Si  0.12500  0.12500  0.12500  1.0  0.006
""")
    log("Wrote Si CIF (Fd-3m, full symmetry ops)")

    # Write initial .instprm
    instprm_path = os.path.join(work_dir, 'initial.instprm')
    with open(instprm_path, 'w') as f:
        f.write(f"""#GSAS-II instrument parameter file; do not add/delete items!
Type:PXC
Lam1:{LAM1}
Lam2:{LAM2}
I(L2)/I(L1):{I_RATIO}
Zero:0.0
Polariz.:{POLARIZ}
U:2.0
V:-2.0
W:5.0
X:0.5
Y:0.5
Z:0.0
SH/L:0.01
Azimuth:0.0
""")
    log("Wrote initial .instprm")

    # ── Import GSAS-II ─────────────────────────────────────────────────
    log("Importing GSAS-II...")
    _add_gsas2_paths()
    try:
        import GSASIIscriptable as G2sc
    except ImportError:
        try:
            from GSASII import GSASIIscriptable as G2sc
        except ImportError:
            log("ERROR: Cannot import GSASIIscriptable")
            sys.exit(1)
    log("GSAS-II imported.")

    # ── Create project ─────────────────────────────────────────────────
    gpx_path = os.path.join(work_dir, 'si_cal.gpx')
    log("Creating GSAS-II project...")
    gpx = G2sc.G2Project(newgpx=gpx_path)
    log("Project created.")

    # Add histogram
    log("Adding histogram...")
    hist = gpx.add_powder_histogram(xye_path, instprm_path,
                                     fmthint='xye')
    if isinstance(hist, list):
        hist = hist[0]
    log(f"Histogram added.")

    # Fix histogram scale
    try:
        hist.data['Sample Parameters']['Scale'] = [1.0, False]
    except Exception as e:
        log(f"  Warning: could not fix hist scale: {e}")

    # Add phase
    log("Adding Si phase from CIF...")
    try:
        phase = gpx.add_phase(cif_path, phasename='Si_640g',
                               histograms=[hist])
        if isinstance(phase, list):
            phase = phase[0]
        log("Phase added successfully.")
    except Exception as e:
        log(f"ERROR adding phase: {e}")
        traceback.print_exc()
        sys.exit(1)

    # Fix cell (NIST certified, don't refine)
    phase.set_refinements({'Cell': False})

    # Set initial scale
    try:
        hapData = list(phase.data['Histograms'].values())[0]
        hapData['Scale'] = [1.0, True]
    except Exception as e:
        log(f"  Warning: could not set phase scale: {e}")

    gpx.save()
    log("Project saved. Starting refinement stages...\n")

    # ── Stage 1: Background + scale ───────────────────────────────────
    log("Stage 1: Background + scale...")
    try:
        gpx.do_refinements([{
            'set': {
                'Background': {'type': 'chebyschev-1', 'refine': True,
                                'no. coeffs': 8},
            },
            'cycles': 5,
        }])
        stats = hist.get_statistics()
        log(f"  Rwp = {stats.get('Rwp', '?')}%")
    except Exception as e:
        log(f"  Stage 1 failed: {e}")
        traceback.print_exc()

    # ── Stage 2: Profile (U, V, W, X, Y, Zero, SH/L) ─────────────────
    log("\nStage 2: Profile + Zero + SH/L...")
    try:
        gpx.do_refinements([{
            'set': {
                'Background': {'type': 'chebyschev-1', 'refine': True,
                                'no. coeffs': 8},
                'Instrument Parameters': ['U', 'V', 'W', 'X', 'Y',
                                           'Zero', 'SH/L'],
            },
            'cycles': 15,
        }])
        stats = hist.get_statistics()
        log(f"  Rwp = {stats.get('Rwp', '?')}%")
    except Exception as e:
        log(f"  Stage 2 failed: {e}")
        traceback.print_exc()

    # ── Stage 3: Final polish ──────────────────────────────────────────
    log("\nStage 3: Final polish (20 cycles)...")
    try:
        gpx.do_refinements([{
            'set': {
                'Background': {'type': 'chebyschev-1', 'refine': True,
                                'no. coeffs': 8},
                'Instrument Parameters': ['U', 'V', 'W', 'X', 'Y',
                                           'Zero', 'SH/L'],
            },
            'cycles': 20,
        }])
        stats = hist.get_statistics()
        log(f"  Rwp = {stats.get('Rwp', '?')}%")
    except Exception as e:
        log(f"  Stage 3 failed: {e}")
        traceback.print_exc()

    # ── Extract refined parameters ─────────────────────────────────────
    log("\nExtracting refined instrument parameters...")
    inst = hist.data['Instrument Parameters'][0]

    params = {}
    for key in ['U', 'V', 'W', 'X', 'Y', 'Zero', 'SH/L']:
        entry = inst.get(key, [0, 0])
        params[key] = float(entry[1]) if len(entry) > 1 else float(entry[0])

    log("\n" + "=" * 60)
    log("REFINED INSTRUMENT PARAMETERS (GSAS-II TCH profile)")
    log("=" * 60)
    for k, v in params.items():
        log(f"  {k:6s} = {v:.6f}")

    # ── Write .instprm ────────────────────────────────────────────────
    instprm_content = f"""#GSAS-II instrument parameter file; do not add/delete items!
Type:PXC
Lam1:{LAM1}
Lam2:{LAM2}
I(L2)/I(L1):{I_RATIO}
Zero:{params['Zero']:.5f}
Polariz.:{POLARIZ}
U:{params['U']:.4f}
V:{params['V']:.4f}
W:{params['W']:.4f}
X:{params['X']:.4f}
Y:{params['Y']:.4f}
Z:0.0
SH/L:{params['SH/L']:.5f}
Azimuth:0.0
"""

    with open(OUTPUT_INSTPRM, 'w') as f:
        f.write(instprm_content)

    log(f"\nWritten to: {OUTPUT_INSTPRM}")
    log(instprm_content)

    # ── Cleanup ────────────────────────────────────────────────────────
    try:
        shutil.rmtree(work_dir)
    except Exception:
        log(f"(temp dir not cleaned: {work_dir})")

    log("Done!")


def _add_gsas2_paths():
    prefix = sys.prefix
    candidates = [
        os.path.join(prefix, 'GSAS-II', 'GSASII'),
        os.path.join(prefix, 'GSAS-II'),
        os.path.join(prefix, 'GSAS-II', 'backcompat'),
    ]
    for p in candidates:
        if os.path.isdir(p) and p not in sys.path:
            sys.path.insert(0, p)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f"\nFATAL ERROR: {e}", flush=True)
        traceback.print_exc()
        sys.exit(1)
