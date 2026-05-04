"""Si 640g calibration — v5 (no add_phase, direct data injection)"""
import os, sys, tempfile, shutil, traceback, copy
import numpy as np
print("=== Si 640g Calibration v5 ===", flush=True)

HERE = os.path.dirname(os.path.abspath(__file__))
XY = os.path.join(HERE, 'POROSOFF_NIST-Si640',
    'BB_0.01st-12sp-0.5IS-20RS1-OpenRS2--NIST-Si640_001.xy')
OUT = os.path.join(HERE, 'smartlab_Si640g.instprm')

if not os.path.isfile(XY):
    print(f"ERROR: {XY} not found", flush=True); sys.exit(1)

data = np.loadtxt(XY)
tt, intensity = data[:,0], data[:,1]
sigma = np.sqrt(np.maximum(intensity, 1.0))
print(f"Data: {len(tt)} pts, {tt[0]:.1f}-{tt[-1]:.1f} deg", flush=True)

W = tempfile.mkdtemp(prefix='si_cal_')
xye = os.path.join(W, 'si.xye')
np.savetxt(xye, np.column_stack([tt, intensity, sigma]), fmt='%.6f  %.4f  %.4f')

ins = os.path.join(W, 'start.instprm')
with open(ins, 'w') as f:
    f.write("""#GSAS-II instrument parameter file; do not add/delete items!
Type:PXC
Lam1:1.540593
Lam2:1.544414
I(L2)/I(L1):0.5
Zero:0.0
Polariz.:0.5
U:2.0
V:-2.0
W:5.0
X:0.5
Y:0.5
Z:0.0
SH/L:0.01
Azimuth:0.0
""")
print("Setup done", flush=True)

# Import GSAS-II
for p in [os.path.join(sys.prefix, 'GSAS-II', 'GSASII'),
          os.path.join(sys.prefix, 'GSAS-II'),
          os.path.join(sys.prefix, 'GSAS-II', 'backcompat')]:
    if os.path.isdir(p) and p not in sys.path:
        sys.path.insert(0, p)
try:
    import GSASIIscriptable as G2sc
except ImportError:
    from GSASII import GSASIIscriptable as G2sc
print("GSAS-II ready", flush=True)

# Also import internal modules we need
try:
    import GSASIIobj as G2obj
    import GSASIIspc as G2spc
    import GSASIIlattice as G2lat
    import GSASIIElem as G2elem
except ImportError:
    from GSASII import GSASIIobj as G2obj
    from GSASII import GSASIIspc as G2spc
    from GSASII import GSASIIlattice as G2lat
    from GSASII import GSASIIElem as G2elem
print("GSAS-II internals loaded", flush=True)

# Project + histogram
gpx_path = os.path.join(W, 'cal.gpx')
gpx = G2sc.G2Project(newgpx=gpx_path)
hist = gpx.add_powder_histogram(xye, ins, fmthint='xye')
if isinstance(hist, list): hist = hist[0]
histname = hist.name
print(f"Histogram: {histname}", flush=True)

# ── Build Si phase by direct injection into gpx.data ──────────────────
print("Building Si phase directly (bypassing add_phase)...", flush=True)

# 1. Get space group data
print("  Loading space group F m -3 m (#225)...", flush=True)
sgErr, SGData = G2spc.SpcGroup('F m -3 m')
if sgErr:
    print(f"  SpcGroup error: {sgErr}", flush=True)
    sys.exit(1)
print(f"  SG: {SGData.get('SpGrp','?')}", flush=True)

# 2. Build cell
a = 5.431109
cellList = [False, a, a, a, 90.0, 90.0, 90.0, a**3]  # [refine, a,b,c,alpha,beta,gamma,vol]
print(f"  Cell: a={a} A, V={a**3:.2f} A^3", flush=True)

# 3. Create phase using SetNewPhase
print("  Creating phase data structure...", flush=True)
phasename = 'Si_cal'
phaseData = G2obj.SetNewPhase(Name=phasename, SGData=SGData, cell=cellList)
print("  Phase data created", flush=True)

# 4. Add atoms: Si at (0,0,0) and (1/4,1/4,1/4) for diamond structure
#    Using Fm-3m with two sites to emulate Fd-3m reflections
print("  Adding Si atoms...", flush=True)
generalData = phaseData['General']
generalData['AtomTypes'] = ['Si']
generalData['NoAtoms'] = {'Si': 2}
generalData['Type'] = 'nuclear'

# Get atom info for Si
atomInfo = G2elem.GetAtomInfo('Si')
print(f"  Si atom info loaded", flush=True)

# Build atom entries. GSAS-II atom list format:
# [label, type, refflags, x, y, z, frac, Utype, Uiso, Id, ...]
# The exact format is version-dependent; let's build minimal entries
atomData = phaseData['Atoms']

for i, (lbl, x, y, z) in enumerate([('Si1', 0.0, 0.0, 0.0),
                                      ('Si2', 0.25, 0.25, 0.25)]):
    # Use MakeNewAtom if available, else build manually
    try:
        atom = [lbl, 'Si', '', x, y, z, 1.0, 'Uiso', 0.006, i]
        atomData.append(atom)
    except Exception as e:
        print(f"  Atom add error: {e}", flush=True)

print(f"  {len(atomData)} atoms added", flush=True)

# 5. Create histogram-phase (HAP) linkage
print("  Creating histogram-phase link...", flush=True)
hapKey = histname
phaseData['Histograms'][hapKey] = {
    'Babinet': {'BabA': [0.0, False], 'BabU': [0.0, False]},
    'Extinction': [0.0, False],
    'Flack': [0.0, False],
    'HStrain': [np.zeros(6, dtype=float).tolist(), np.zeros(6, dtype=bool).tolist()],
    'Mustrain': ['isotropic', [1000.0, False], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]],
    'Pref.Ori.': ['MD', 1.0, False, [0, 0, 1], 0, [], {}],
    'Scale': [1.0, True],
    'Size': ['isotropic', [1.0, False], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]],
    'Use': True,
    'Fix FXU': ' ',
}
print("  HAP link created", flush=True)

# 6. Inject phase into project
print("  Injecting phase into project...", flush=True)

# Find where phases are stored in gpx.data
# GSAS-II stores phases in gpx.data with key 'Phases'
if 'Phases' not in gpx.data:
    gpx.data['Phases'] = {}
gpx.data['Phases'][phasename] = phaseData

# Also need to register in the histogram's phase list
try:
    histData = hist.data
    if 'Reflection Lists' not in histData:
        histData['Reflection Lists'] = {}
    histData['Reflection Lists'][phasename] = {}
except Exception as e:
    print(f"  RefList setup warning: {e}", flush=True)

gpx.save()
print("  Phase injected and saved", flush=True)

# Verify phase is visible
phases = gpx.phases()
print(f"  Project has {len(phases)} phase(s)", flush=True)
if len(phases) == 0:
    print("ERROR: Phase not visible to GSAS-II after injection", flush=True)
    print("  gpx.data keys:", list(gpx.data.keys()), flush=True)
    sys.exit(1)

phase = phases[0]
phase.set_refinements({'Cell': False})
print("Phase ready\n", flush=True)

# Refinement
print("--- Stage 1: BG + scale ---", flush=True)
try:
    gpx.do_refinements([{'set': {
        'Background': {'type': 'chebyschev-1', 'refine': True, 'no. coeffs': 8},
    }, 'cycles': 5}])
    print(f"  Rwp = {hist.get_statistics().get('Rwp','?')}%", flush=True)
except Exception as e:
    print(f"  Stage 1 error: {e}", flush=True)
    traceback.print_exc()

print("\n--- Stage 2: Profile ---", flush=True)
try:
    gpx.do_refinements([{'set': {
        'Background': {'type': 'chebyschev-1', 'refine': True, 'no. coeffs': 8},
        'Instrument Parameters': ['U','V','W','X','Y','Zero','SH/L'],
    }, 'cycles': 15}])
    print(f"  Rwp = {hist.get_statistics().get('Rwp','?')}%", flush=True)
except Exception as e:
    print(f"  Stage 2 error: {e}", flush=True)
    traceback.print_exc()

print("\n--- Stage 3: Polish ---", flush=True)
try:
    gpx.do_refinements([{'set': {
        'Background': {'type': 'chebyschev-1', 'refine': True, 'no. coeffs': 8},
        'Instrument Parameters': ['U','V','W','X','Y','Zero','SH/L'],
    }, 'cycles': 20}])
    print(f"  Rwp = {hist.get_statistics().get('Rwp','?')}%", flush=True)
except Exception as e:
    print(f"  Stage 3 error: {e}", flush=True)
    traceback.print_exc()

# Extract
inst = hist.data['Instrument Parameters'][0]
params = {}
for k in ['U','V','W','X','Y','Zero','SH/L']:
    e = inst.get(k, [0,0])
    params[k] = float(e[1]) if len(e) > 1 else float(e[0])

print("\n" + "="*50, flush=True)
print("REFINED PARAMETERS", flush=True)
print("="*50, flush=True)
for k,v in params.items():
    print(f"  {k:6s} = {v:.6f}", flush=True)

content = f"""#GSAS-II instrument parameter file; do not add/delete items!
Type:PXC
Lam1:1.540593
Lam2:1.544414
I(L2)/I(L1):0.5
Zero:{params['Zero']:.5f}
Polariz.:0.5
U:{params['U']:.4f}
V:{params['V']:.4f}
W:{params['W']:.4f}
X:{params['X']:.4f}
Y:{params['Y']:.4f}
Z:0.0
SH/L:{params['SH/L']:.5f}
Azimuth:0.0
"""
with open(OUT, 'w') as f:
    f.write(content)
print(f"\nSaved: {OUT}", flush=True)
print(content, flush=True)
try: shutil.rmtree(W)
except: pass
print("DONE", flush=True)
