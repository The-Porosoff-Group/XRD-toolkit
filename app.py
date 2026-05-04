"""
app.py  —  Catalysis Data Toolkit
Local web server. Run with:  python app.py
Then open:  http://localhost:5000
"""

import os, sys, re, json, base64, webbrowser

# Force unbuffered/line-buffered stdout so GSAS-II refinement progress
# appears in terminal immediately.  os.environ alone doesn't work because
# Python's IO is already initialized by the time app.py runs.
os.environ['PYTHONUNBUFFERED'] = '1'
try:
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)
except AttributeError:
    pass  # Python < 3.7
from datetime import datetime
from threading import Timer

# ── Dependency check ─────────────────────────────────────────────────────────
_REQUIRED = {
    'numpy':      'numpy',
    'yaml':       'pyyaml',
    'flask':      'flask',
    'scipy':      'scipy',
    'matplotlib': 'matplotlib',
    'requests':   'requests',
    'pandas':     'pandas',
}
_missing = []
for _mod, _pkg in _REQUIRED.items():
    try:
        __import__(_mod)
    except ImportError:
        _missing.append(_pkg)

if _missing:
    print()
    print("  ╔══════════════════════════════════════════════════════╗")
    print("  ║  Missing required packages:                         ║")
    for _pkg in _missing:
        print(f"  ║    - {_pkg:<49}║")
    print("  ║                                                     ║")
    print("  ║  Option 1 (recommended): Use run.bat instead        ║")
    print("  ║    It sets up everything automatically.              ║")
    print("  ║                                                     ║")
    print("  ║  Option 2: Activate the conda env first, then run:  ║")
    print("  ║    conda activate .conda_env                        ║")
    print("  ║    python app.py                                    ║")
    print("  ║                                                     ║")
    print("  ║  Option 3: Install manually:                        ║")
    print(f"  ║    pip install {' '.join(_missing):<39}║")
    print("  ╚══════════════════════════════════════════════════════╝")
    print()
    input("  Press Enter to exit...")
    sys.exit(1)

import numpy as np
import yaml
from flask import Flask, render_template, request, jsonify, send_file
from flask.json.provider import DefaultJSONProvider


class NumpyJSONProvider(DefaultJSONProvider):
    """JSON provider that serialises numpy scalars and arrays."""

    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)

# ── Path setup ────────────────────────────────────────────────────────────────
BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
MODULES_DIR = os.path.join(BASE_DIR, 'modules')
UPLOAD_DIR  = os.path.join(BASE_DIR, 'uploads')
CONFIG_PATH = os.path.join(BASE_DIR, 'config.yaml')
XRD_PRESETS_PATH = os.path.join(BASE_DIR, 'xrd_refinement_presets.json')
os.makedirs(UPLOAD_DIR, exist_ok=True)
sys.path.insert(0, MODULES_DIR)


def _normalize_formula_case(formula):
    """Allow all-lowercase formula input for formula searches."""
    formula = (formula or '').strip().replace(' ', '')
    if not formula or re.search(r'[A-Z]', formula):
        return formula
    valid = {
        'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si',
        'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',
        'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb',
        'Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
        'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho',
        'Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg',
        'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np',
        'Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg',
        'Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og',
    }
    out = []
    i = 0
    while i < len(formula):
        if formula[i].isdigit():
            out.append(formula[i])
            i += 1
            continue
        if not formula[i].isalpha():
            out.append(formula[i])
            i += 1
            continue
        two = formula[i:i+2].capitalize()
        one = formula[i].upper()
        if i + 1 < len(formula) and two in valid:
            out.append(two)
            i += 2
        elif one in valid:
            out.append(one)
            i += 1
        else:
            out.append(formula[i].upper())
            i += 1
    return ''.join(out)

# ── Load config ───────────────────────────────────────────────────────────────
def load_config():
    defaults = {
        'materials_project': {'api_key': ''},
        'cache': {'directory': '~/.catalysis_toolkit_cache', 'max_size_mb': 500},
        'performance': {'max_outer_iterations': 10, 'preload_pymatgen': True},
    }
    if os.path.exists(CONFIG_PATH):
        try:
            with open(CONFIG_PATH) as f:
                user = yaml.safe_load(f) or {}
            # Deep merge
            for section, vals in user.items():
                if isinstance(vals, dict):
                    defaults.setdefault(section, {}).update(vals)
                else:
                    defaults[section] = vals
        except Exception as e:
            print(f"  Warning: could not read config.yaml: {e}")
    return defaults

CONFIG = load_config()
MP_API_KEY   = CONFIG['materials_project'].get('api_key', '').strip()
CACHE_DIR    = CONFIG['cache'].get('directory', '~/.catalysis_toolkit_cache')
CACHE_MAX_MB = CONFIG['cache'].get('max_size_mb', 500)
MAX_OUTER    = CONFIG['performance'].get('max_outer_iterations', 10)

# ── Import modules ────────────────────────────────────────────────────────────
import gc_processor
import tga_processor
import bet_processor
import xrd_processor
from xrd.cif_cache import get_cache, cached_fetch_cod, cached_fetch_mp
from xrd.mp_api    import (search_by_elements  as mp_search_elements,
                            search_by_formula   as mp_search_formula,
                            search_by_name      as mp_search_name,
                            fetch_cif           as mp_fetch_cif,
                            validate_api_key    as mp_validate_key)
from xrd.cod_api   import (search_by_elements  as cod_search_elements,
                            search_by_formula   as cod_search_formula,
                            search_by_name      as cod_search_name,
                            fetch_cif           as cod_fetch_cif,
                            get_stick_pattern,  SORT_OPTIONS)

# Initialise cache with config settings
_cache = get_cache(cache_dir=CACHE_DIR, max_size_mb=CACHE_MAX_MB)

# ── Preload pymatgen ──────────────────────────────────────────────────────────
_pymatgen_ready = False
if CONFIG['performance'].get('preload_pymatgen', True):
    try:
        print("  Pre-loading pymatgen...", end=' ', flush=True)
        from pymatgen.core import Structure
        from pymatgen.analysis.diffraction.xrd import XRDCalculator
        _dummy = XRDCalculator(wavelength=1.54056)
        _pymatgen_ready = True
        print("done.")
    except ImportError:
        print("not installed (run: pip install pymatgen)")
    except Exception as e:
        print(f"failed ({e})")

# ── Flask app ─────────────────────────────────────────────────────────────────
app = Flask(__name__)
app.json_provider_class = NumpyJSONProvider
app.json = NumpyJSONProvider(app)
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024

MODULES = [
    {'id': 'gc',  'name': 'GC Analysis',       'icon': '⚗️',  'status': 'active'},
    {'id': 'tga', 'name': tga_processor.MODULE_INFO['name'],
                  'icon': tga_processor.MODULE_INFO['icon'],
                  'status': tga_processor.MODULE_INFO['status']},
    {'id': 'bet', 'name': bet_processor.MODULE_INFO['name'],
                  'icon': bet_processor.MODULE_INFO['icon'],
                  'status': bet_processor.MODULE_INFO['status']},
    {'id': 'xrd', 'name': xrd_processor.MODULE_INFO['name'],
                  'icon': xrd_processor.MODULE_INFO['icon'],
                  'status': 'active'},
]

# ── Routes — General ──────────────────────────────────────────────────────────

@app.route('/')
def index():
    reaction_configs = gc_processor.list_reaction_configs(
        os.path.join(MODULES_DIR, 'reaction_configs'))
    return render_template('xrd_toolkit/index.html',
                           modules=MODULES,
                           reaction_configs=reaction_configs,
                           mp_key_set=bool(MP_API_KEY),
                           pymatgen_ready=_pymatgen_ready)


@app.route('/xrd')
@app.route('/xrd-toolkit')
def xrd_toolkit():
    reaction_configs = gc_processor.list_reaction_configs(
        os.path.join(MODULES_DIR, 'reaction_configs'))
    return render_template('xrd_toolkit/index.html',
                           modules=MODULES,
                           reaction_configs=reaction_configs,
                           mp_key_set=bool(MP_API_KEY),
                           pymatgen_ready=_pymatgen_ready)


@app.route('/api/reaction_configs')
def get_reaction_configs():
    return jsonify(gc_processor.list_reaction_configs(
        os.path.join(MODULES_DIR, 'reaction_configs')))

@app.route('/api/status')
def api_status():
    cache_stats = _cache.stats()
    return jsonify({
        'pymatgen_ready': _pymatgen_ready,
        'mp_key_set':     bool(MP_API_KEY),
        'cache':          cache_stats,
    })

@app.route('/api/cache/clear', methods=['POST'])
def clear_cache():
    _cache.clear()
    return jsonify({'ok': True, 'message': 'Cache cleared.'})

# ── Routes — GC ───────────────────────────────────────────────────────────────

@app.route('/api/process_gc', methods=['POST'])
def process_gc():
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'No file uploaded.'}), 400
        f = request.files['file']
        if not f.filename.endswith('.xlsx'):
            return jsonify({'error': 'Please upload an .xlsx file.'}), 400
        safe_name   = re.sub(r'[^\w\-.]', '_', f.filename)
        upload_path = os.path.join(UPLOAD_DIR, safe_name)
        f.save(upload_path)
        form        = request.form
        config_file = form.get('reaction_config')
        config_path = os.path.join(MODULES_DIR, 'reaction_configs', config_file)
        if not os.path.exists(config_path):
            return jsonify({'error': f'Reaction config not found: {config_file}'}), 400
        reaction_config = gc_processor.load_reaction_config(config_path)
        metadata = {
            'catalyst_id': form.get('catalyst_id', 'Unknown'),
            'reactant':    reaction_config['reactant'],
            'temperature': form.get('temperature', ''),
            'pressure':    form.get('pressure', ''),
            'ghsv':        form.get('ghsv', ''),
            'notes':       form.get('notes', ''),
            'source_file': f.filename,
            'reaction':    reaction_config['name'],
        }
        inlet_flows = {}
        try:
            inlet_flows = {k: float(v)
                           for k, v in json.loads(form.get('inlet_flows', '{}')).items()
                           if v}
        except (json.JSONDecodeError, ValueError):
            pass
        ss_start = int(form.get('ss_start', 1))
        ss_end   = int(form.get('ss_end',   999))
        output_base = form.get('output_dir', '').strip()
        if not output_base or not os.path.isdir(output_base):
            output_base = os.path.join(BASE_DIR, 'results')
        ts      = datetime.now().strftime('%Y%m%d_%H%M%S')
        safe_id = re.sub(r'[^\w\-]', '_', metadata['catalyst_id'])
        output_dir = os.path.join(output_base, f'{safe_id}_{ts}')
        result = gc_processor.run(
            filepath=upload_path, output_dir=output_dir,
            reaction_config=reaction_config, metadata=metadata,
            inlet_flows=inlet_flows, ss_start=ss_start, ss_end=ss_end)
        with open(result['plot_path'], 'rb') as img:
            plot_b64 = base64.b64encode(img.read()).decode('utf-8')
        clean = {k: v for k, v in result.items() if k != 'plot_path'}
        clean['plot_b64'] = plot_b64
        return jsonify(clean)
    except Exception as e:
        import traceback
        return jsonify({'error': str(e), 'trace': traceback.format_exc()}), 500

@app.route('/api/download')
def download_file():
    path = request.args.get('path', '')
    if os.path.isfile(path):
        return send_file(path, as_attachment=True)
    return jsonify({'error': 'File not found'}), 404

@app.route('/api/open_folder')
def open_folder():
    path = request.args.get('path', '')
    if os.path.isdir(path):
        os.startfile(path)
        return jsonify({'ok': True})
    return jsonify({'error': 'Folder not found'}), 404

# ── Routes — XRD search ───────────────────────────────────────────────────────

@app.route('/api/xrd/search', methods=['POST'])
def xrd_search():
    try:
        data       = request.get_json()
        elements   = [e.strip() for e in data.get('elements', []) if e.strip()]
        name       = data.get('name', '').strip()
        formula    = _normalize_formula_case(data.get('formula', ''))
        wavelength = float(data.get('wavelength', 1.54056))
        sort_by    = data.get('sort_by', 'formula')
        strict     = bool(data.get('strict', True))
        limit      = min(int(data.get('limit', 100)), 200)
        source     = data.get('source', 'cod')  # 'cod' or 'mp' or 'both'

        cod_results = []
        mp_results  = []

        # ── COD search ────────────────────────────────────────────────────────
        if source in ('cod', 'both'):
            if formula:
                r = cod_search_formula(formula, max_results=limit, sort_by=sort_by)
            elif name:
                r = cod_search_name(name, max_results=limit, sort_by=sort_by)
            elif elements:
                r = cod_search_elements(elements, strict=strict,
                                         max_results=limit, sort_by=sort_by)
            else:
                r = {'error': 'Provide elements, formula, or name.'}
            if isinstance(r, list):
                for entry in r:
                    entry['source'] = 'cod'
                cod_results = r
            elif isinstance(r, dict) and 'error' in r and source == 'cod':
                return jsonify(r), 500

        # ── Materials Project search ──────────────────────────────────────────
        if source in ('mp', 'both') and MP_API_KEY:
            if formula:
                r = mp_search_formula(formula, MP_API_KEY,
                                       max_results=limit, sort_by=sort_by)
            elif elements:
                r = mp_search_elements(elements, MP_API_KEY, strict=strict,
                                        max_results=limit, sort_by=sort_by)
            elif name:
                # MP has no free-text search — use our smart name→elements fallback
                r = mp_search_name(name, MP_API_KEY,
                                    max_results=limit, sort_by=sort_by)
            else:
                r = []
            if isinstance(r, list):
                mp_results = r
            elif isinstance(r, dict) and 'error' in r and source == 'mp':
                return jsonify(r), 500

        # Combine: MP first (computed, more complete), then COD
        combined = mp_results + cod_results

        # Add stick patterns (first 20 only — performance)
        # Validate each entry first so SG number and system are correct
        try:
            from xrd import validate_phases as _vp
        except ImportError:
            _vp = None

        for entry in combined[:20]:
            try:
                ph_seed = dict(entry)
                # MP search results carry inline CIF as _cif_text;
                # get_stick_pattern expects cif_text.
                if ph_seed.get('_cif_text') and not ph_seed.get('cif_text'):
                    ph_seed['cif_text'] = ph_seed['_cif_text']
                if _vp:
                    validated = _vp([dict(ph_seed)], fetch_missing=False)
                    ph = validated[0] if validated else ph_seed
                else:
                    ph = ph_seed
                # Preserve CIF if validate_phases dropped it
                if ph_seed.get('cif_text') and not ph.get('cif_text'):
                    ph['cif_text'] = ph_seed['cif_text']
                entry['stick_pattern'] = get_stick_pattern(ph, wavelength)
            except Exception as _e:
                print(f"  Stick pattern failed for "
                      f"{entry.get('formula', '?')}: {_e}", flush=True)
                entry['stick_pattern'] = []

        # Cache any CIF text that came back inline with MP results
        for entry in combined:
            cif_text = entry.pop('_cif_text', '')
            if cif_text and '_cell_length_a' in cif_text:
                mp_id     = entry.get('mp_id', entry.get('cod_id', ''))
                cache_key = f"mp:{mp_id}"
                if mp_id and not _cache.has(cache_key):
                    _cache.put(cache_key, cif_text)

        return jsonify({
            'results':     combined,
            'total':       len(combined),
            'cod_count':   len(cod_results),
            'mp_count':    len(mp_results),
            'mp_available': bool(MP_API_KEY),
        })

    except Exception as e:
        import traceback
        return jsonify({'error': str(e), 'trace': traceback.format_exc()}), 500


@app.route('/api/xrd/fetch_cif', methods=['POST'])
def xrd_fetch_cif():
    """Fetch and cache a CIF from COD or Materials Project."""
    try:
        data       = request.get_json()
        cod_id     = str(data.get('cod_id', ''))
        source_raw = data.get('source', 'cod')
        source     = ('mp' if source_raw in ('mp', 'materials_project')
                      or cod_id.startswith('mp-') else 'cod')
        wavelength = float(data.get('wavelength', 1.54056))

        phase_hint = data.get('phase_hint') or {}

        if source == 'mp':
            result = cached_fetch_mp(cod_id, MP_API_KEY, mp_fetch_cif)
        else:
            result = cached_fetch_cod(cod_id, cod_fetch_cif)

        # ── Diagnostic: confirm W2C mp-2034 is using the Pbcn fixture ──
        if str(cod_id) == 'mp-2034':
            _ct = result.get('cif_text', '') or ''
            _p1_marker = "_symmetry_space_group_name_H-M   'P 1'"
            _has_pbcn = ('Pbcn' in _ct)
            _has_p1   = (_p1_marker in _ct) or ("'P 1'" in _ct and not _has_pbcn)
            print(f"  W2C mp-2034 CIF source check: "
                  f"SG={result.get('spacegroup_number')} "
                  f"formula={result.get('formula')} "
                  f"has_Pbcn={_has_pbcn} has_P1={_has_p1}",
                  flush=True)

        # Merge search-row metadata into fetched result.
        # MP fetch often returns P1/full-cell (SG=1, formula=W8C4)
        # while the search correctly had SG=60, formula=W2C.
        #
        # CRITICAL: when a fixture is providing the CIF (e.g. mp-2034
        # W2C → fixtures/w2c_pbcn_mp_2034.cif), the fixture's cell is in
        # the CifWriter / Pbcn convention, but phase_hint['a','b','c']
        # comes from parse_cif's sorted-axis convention.  Overriding
        # the fixture cell with phase_hint axes produces a CIF whose
        # cell and atomic positions are in *different* conventions —
        # exactly the disaster that broke W2C refinement.
        #
        # Allowed phase_hint overrides: labels and intent (formula,
        # name, SG symbol/number, crystal system).
        # Forbidden for MP CIFs: cell axes and angles. Those must come
        # from the CIF that also provides the atom positions.
        try:
            from xrd.mp_api import _fixture_cif_for as _fix_lookup
            _has_fixture = bool(_fix_lookup(str(cod_id)))
        except Exception:
            _has_fixture = False

        if phase_hint and source == 'mp':
            _allowed_keys = ['formula', 'name', 'spacegroup',
                              'spacegroup_number', 'system']
            for _hk in _allowed_keys:
                _hv = phase_hint.get(_hk)
                if _hv not in (None, '', 0):
                    result[_hk] = _hv
            print(f"  phase_hint merge: applied display/intent metadata "
                  f"only for {cod_id}; CIF cell remains authoritative.",
                  flush=True)

        # Store CIF text server-side; don't send over wire
        cif_text  = result.pop('cif_text', '')
        cache_key = f"{'mp' if source=='mp' else 'cod'}:{cod_id}"
        _cache.put(cache_key, cif_text)
        if source == 'mp':
            _cache.put(f"{cache_key}:gsas:v1", cif_text)

        cif_check = {
            'status': 'error' if not cif_text else 'ok',
            'ok': bool(cif_text),
            'messages': [],
            'expected_sg': result.get('spacegroup_number'),
            'parsed_sg': None,
            'site_count': None,
            'cell': None,
        }
        if cif_text:
            try:
                from modules.xrd.crystallography import parse_cif as _parse_cif
                _parsed_cif = _parse_cif(cif_text)
                _parsed_sg = int(_parsed_cif.get('spacegroup_number', 1) or 1)
                _expected_sg = int(result.get('spacegroup_number', 0) or 0)
                _sites = _parsed_cif.get('sites') or []
                cif_check.update({
                    'parsed_sg': _parsed_sg,
                    'site_count': len(_sites),
                    'cell': {
                        'a': _parsed_cif.get('a'),
                        'b': _parsed_cif.get('b'),
                        'c': _parsed_cif.get('c'),
                        'alpha': _parsed_cif.get('alpha'),
                        'beta': _parsed_cif.get('beta'),
                        'gamma': _parsed_cif.get('gamma'),
                    },
                })
                if _expected_sg > 1 and _parsed_sg == 1:
                    cif_check['status'] = 'warn'
                    cif_check['ok'] = True
                    cif_check['messages'].append(
                        'Source CIF is P1/full-cell; GSAS-II prep will '
                        'attempt a symmetry-safe reduction before refinement.')
                elif (_expected_sg > 1 and _parsed_sg > 1
                        and _parsed_sg != _expected_sg):
                    cif_check['status'] = 'error'
                    cif_check['ok'] = False
                    cif_check['messages'].append(
                        f"Source CIF declares SG {_parsed_sg}, but the "
                        f"selected phase expects SG {_expected_sg}.")
                elif _sites:
                    cif_check['messages'].append(
                        f"CIF parsed as SG {_parsed_sg} with "
                        f"{len(_sites)} atom site(s).")
                else:
                    cif_check['status'] = 'warn'
                    cif_check['ok'] = True
                    cif_check['messages'].append(
                        'CIF parsed, but no atom-site loop was found.')
            except Exception as _ve:
                cif_check['status'] = 'error'
                cif_check['ok'] = False
                cif_check['messages'].append(f'CIF parse check failed: {_ve}')

        # Compute accurate stick pattern server-side using crystallography engine
        # IMPORTANT: keep CIF text available for stick generation —
        # it was popped from result above, so re-attach it.
        # Also preserve the original SG/formula from the MP search result,
        # because validate_phases may re-parse the P1 CIF and overwrite
        # SG to 1 and formula to the full-cell formula (e.g. W8C4).
        try:
            from xrd import validate_phases
            ph_seed = dict(result)
            ph_seed['cif_text'] = cif_text
            ph_seed['source'] = source
            validated = validate_phases([dict(ph_seed)], fetch_missing=False)
            ph_for_sticks = validated[0] if validated else ph_seed
            # Preserve CIF if validate_phases dropped it
            if cif_text and not ph_for_sticks.get('cif_text'):
                ph_for_sticks['cif_text'] = cif_text
            # Restore search-row metadata if validation downgraded to P1
            if (ph_seed.get('spacegroup_number', 1) > 1
                    and ph_for_sticks.get('spacegroup_number', 1) == 1):
                for _rk in ('formula', 'spacegroup_number', 'spacegroup',
                             'system', 'source'):
                    if ph_seed.get(_rk):
                        ph_for_sticks[_rk] = ph_seed[_rk]
            sticks = get_stick_pattern(ph_for_sticks, wavelength, tt_min=5.0, tt_max=100.0)
            result['stick_pattern'] = sticks
            result['stick_source'] = 'python_reflections'
        except Exception as _e:
            print(f"  fetch_cif stick pattern failed: {_e}", flush=True)
            result['stick_pattern'] = []
            result['stick_source'] = 'unavailable'

        result['cached'] = result.get('cached', False)
        result['cache_key'] = cache_key
        result['cif_check'] = cif_check
        result['validation'] = cif_check
        result['preview'] = {
            'stick_pattern': result.get('stick_pattern', []),
            'stick_source': result.get('stick_source', 'python_reflections'),
        }
        result['phase'] = {
            k: result.get(k) for k in (
                'formula', 'name', 'spacegroup', 'spacegroup_number',
                'system', 'a', 'b', 'c', 'alpha', 'beta', 'gamma',
                'Z', 'source', 'cod_id', 'mp_id')
            if k in result
        }
        # Ensure source is in the response so the frontend preserves it
        if 'source' not in result:
            result['source'] = source
        return jsonify(result)

    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/xrd/mp_debug_cif', methods=['GET'])
def mp_debug_cif():
    """Show structure data from summary endpoint for a given mp-id."""
    import requests as req
    if not MP_API_KEY:
        return jsonify({'error': 'No MP API key configured'})
    try:
        mp_id   = request.args.get('id', 'mp-91')
        headers = {'X-API-KEY': MP_API_KEY, 'Accept': 'application/json'}
        resp    = req.get(
            'https://api.materialsproject.org/materials/summary/',
            headers=headers,
            params={'material_ids': mp_id, '_fields': 'material_id,formula_pretty,structure',
                    'deprecated': 'false', '_limit': 1},
            timeout=15,
        )
        data = resp.json().get('data', [])
        entry = data[0] if data else {}
        struct = entry.get('structure', {})
        lattice = struct.get('lattice', {}) if struct else {}
        return jsonify({
            'status':          resp.status_code,
            'formula':         entry.get('formula_pretty'),
            'has_structure':   bool(struct),
            'lattice':         lattice,
            'n_sites':         len(struct.get('sites', [])) if struct else 0,
        })
    except Exception as e:
        return jsonify({'error': str(e)})


@app.route('/api/xrd/mp_debug', methods=['GET'])
def mp_debug():
    """Dump raw MP API response for debugging field names."""
    import requests as req
    if not MP_API_KEY:
        return jsonify({'error': 'No MP API key configured'})
    try:
        headers = {'X-API-KEY': MP_API_KEY, 'Accept': 'application/json'}
        resp = req.get('https://api.materialsproject.org/materials/summary/',
                        headers=headers,
                        params={'chemsys': 'W', 'deprecated': 'false', '_limit': 1,
                                '_fields': 'material_id,formula_pretty,symmetry,'
                                           'energy_above_hull,theoretical,structure'},
                        timeout=15)
        return jsonify({
            'status':  resp.status_code,
            'raw':     resp.json(),
        })
    except Exception as e:
        return jsonify({'error': str(e)})


@app.route('/api/xrd/validate_mp_key', methods=['POST'])
def validate_mp_key():
    data = request.get_json() or {}
    key  = data.get('api_key', MP_API_KEY).strip()
    valid, msg = mp_validate_key(key)
    return jsonify({'valid': valid, 'message': msg})


@app.route('/api/xrd/gsas2_status', methods=['GET'])
def gsas2_status():
    """Check if GSAS-II is available in this Python environment."""
    try:
        from modules.xrd.gsasii_backend import is_available, import_error
        return jsonify({
            'available': is_available(),
            'error': import_error(),
        })
    except Exception as e:
        return jsonify({'available': False, 'error': str(e)})


def _load_xrd_presets():
    if not os.path.exists(XRD_PRESETS_PATH):
        return []
    try:
        with open(XRD_PRESETS_PATH, 'r', encoding='utf-8') as f:
            data = json.load(f)
        if isinstance(data, dict):
            presets = data.get('presets', [])
        else:
            presets = data
        return presets if isinstance(presets, list) else []
    except Exception as e:
        print(f"  Warning: could not read XRD presets: {e}", flush=True)
        return []


def _builtin_xrd_presets():
    return [{
        'id': 'builtin-wc-w2c-synergy-s',
        'name': 'WC/W2C Synergy-S production',
        'description': (
            'Built-in recipe for validated WC/W2C GSAS-II fits. Applies '
            'quick constrained fit + cell, phase isolation, zero correction, '
            'free X, Y nonnegative, and WC [001] PO fixed at 0.905.'),
        'builtin': True,
        'locked': True,
        'recipe_only': True,
        'recipe_key': 'wc_w2c_synergy_s',
        'version': 1,
        'phases': [],
        'phase_options': [],
        'controls': {
            'wavelength': '1.54056',
            'wavelength_source': '1.54056',
            'tt_min': '20',
            'tt_max': '60',
            'n_bg_coeffs': 'auto',
            'instrument': 'synergy_s',
            'fix_y_value': '',
            'checkboxes': {
                'xrd-calibration-mode': False,
                'xrd-verification-mode': True,
                'xrd-verify-cell': True,
                'xrd-phase-isolation': True,
                'xrd-zero-not-disp': True,
                'xrd-refine-x': True,
                'xrd-fix-y': False,
                'xrd-y-nonneg': True,
                'xrd-refine-uiso': False,
                'xrd-refine-xyz': False,
            },
        },
        'result_summary': {},
    }]


def _save_xrd_presets(presets):
    tmp_path = XRD_PRESETS_PATH + '.tmp'
    payload = {'schema': 'catalysis-toolkit.xrd-presets.v1',
               'presets': presets}
    with open(tmp_path, 'w', encoding='utf-8') as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    os.replace(tmp_path, XRD_PRESETS_PATH)


def _preset_id_from_name(name):
    slug = re.sub(r'[^a-z0-9]+', '-', (name or '').lower()).strip('-')
    slug = slug[:48] or 'xrd-preset'
    return f"{slug}-{datetime.now().strftime('%Y%m%d%H%M%S')}"


@app.route('/api/xrd/presets', methods=['GET'])
def xrd_list_presets():
    presets = _builtin_xrd_presets() + _load_xrd_presets()
    return jsonify({'presets': presets})


@app.route('/api/xrd/presets', methods=['POST'])
def xrd_save_preset():
    payload = request.get_json(silent=True) or {}
    name = (payload.get('name') or '').strip()
    if not name:
        return jsonify({'error': 'Preset name is required.'}), 400

    presets = _load_xrd_presets()
    preset_id = (payload.get('id') or '').strip()
    if preset_id.startswith('builtin-'):
        return jsonify({'error': 'Built-in presets cannot be overwritten.'}), 400
    now = datetime.now().isoformat(timespec='seconds')

    existing = None
    if preset_id:
        existing = next((p for p in presets if p.get('id') == preset_id), None)
    if existing is None:
        existing = next((p for p in presets if p.get('name') == name), None)

    preset = {
        'id': preset_id or (existing or {}).get('id') or _preset_id_from_name(name),
        'name': name,
        'description': (payload.get('description') or '').strip(),
        'created_at': (existing or {}).get('created_at') or now,
        'updated_at': now,
        'version': 1,
        'phases': payload.get('phases') if isinstance(payload.get('phases'), list) else [],
        'phase_options': (
            payload.get('phase_options')
            if isinstance(payload.get('phase_options'), list) else []),
        'controls': (
            payload.get('controls')
            if isinstance(payload.get('controls'), dict) else {}),
        'result_summary': (
            payload.get('result_summary')
            if isinstance(payload.get('result_summary'), dict) else {}),
    }

    presets = [p for p in presets if p.get('id') != preset['id']]
    presets.append(preset)
    _save_xrd_presets(presets)
    return jsonify({'ok': True, 'preset': preset})


@app.route('/api/xrd/presets/<preset_id>', methods=['DELETE'])
def xrd_delete_preset(preset_id):
    if str(preset_id).startswith('builtin-'):
        return jsonify({'error': 'Built-in presets cannot be deleted.'}), 400
    presets = _load_xrd_presets()
    kept = [p for p in presets if p.get('id') != preset_id]
    if len(kept) == len(presets):
        return jsonify({'error': 'Preset not found.'}), 404
    _save_xrd_presets(kept)
    return jsonify({'ok': True})


@app.route('/api/process_xrd', methods=['POST'])
def process_xrd():
    print("\n=== /api/process_xrd START ===", flush=True)
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'No file uploaded'}), 400
        f = request.files['file']
        safe_name   = re.sub(r'[^\w\-.]', '_', f.filename)
        upload_path = os.path.join(UPLOAD_DIR, safe_name)
        f.save(upload_path)

        form       = request.form
        wavelength = float(form.get('wavelength', 1.54056))
        tt_min     = float(form.get('tt_min', 5.0))
        tt_max     = float(form.get('tt_max', 90.0))
        sample_id  = form.get('sample_id', 'Sample')
        notes      = form.get('notes', '')
        wl_label   = form.get('wavelength_label', f'λ={wavelength:.5f} Å')

        phases = json.loads(form.get('phases', '[]'))
        if not phases:
            return jsonify({'error': 'No phases selected for refinement.'}), 400

        # Re-attach cached CIF text for pymatgen seeding.
        # Also ensure 'source' is preserved — downstream _cif_policy
        # needs it to detect MP phases.
        for ph in phases:
            source = ph.get('source', 'cod')
            ph['source'] = source  # ensure always present
            cid    = str(ph.get('cod_id', ph.get('mp_id', '')))
            # Detect MP phases from ID even if source was lost
            if cid.startswith('mp-') and source == 'cod':
                ph['source'] = 'mp'
                source = 'mp'
            cache_keys = [
                f"{source}:{cid}:gsas:v1",
                f"{source}:{cid}",
                f"cod:{cid}",
                f"mp:{cid}:gsas:v1",
                f"mp:{cid}",
            ]
            text = None
            for cache_key in cache_keys:
                text = _cache.get(cache_key)
                if text:
                    break
            if text:
                ph['cif_text'] = text

        # Optional instrument parameter file for GSAS-II
        instprm_file_path = None
        if 'instprm_file' in request.files:
            instprm_f = request.files['instprm_file']
            if instprm_f.filename:
                safe_instprm = re.sub(r'[^\w\-.]', '_', instprm_f.filename)
                instprm_file_path = os.path.join(UPLOAD_DIR, safe_instprm)
                instprm_f.save(instprm_file_path)

        output_base = form.get('output_dir', '').strip()
        if not output_base or not os.path.isdir(output_base):
            output_base = os.path.join(BASE_DIR, 'results')
        ts      = datetime.now().strftime('%Y%m%d_%H%M%S')
        safe_id = re.sub(r'[^\w\-]', '_', sample_id)
        out_dir = os.path.join(output_base, f'XRD_{safe_id}_{ts}')

        # ── Calibration mode: separate backend ──────────────────────────
        calibration_mode = form.get('calibration_mode', '').lower() == 'true'

        if calibration_mode and form.get('method', '') == 'gsas2':
            from modules.xrd.gsasii_calibration import run_calibration
            from modules.xrd import parse_xrd_file
            from modules.xrd.xrd_plots import make_xrd_plot

            os.makedirs(out_dir, exist_ok=True)

            data = parse_xrd_file(upload_path)
            # Use the first phase as the calibration standard
            if not phases:
                return jsonify({'error': 'No phase selected for calibration.'}), 400
            cal_phase = phases[0]
            # Re-attach CIF text
            source = cal_phase.get('source', 'cod')
            cid = str(cal_phase.get('cod_id', cal_phase.get('mp_id', '')))
            for key_fmt in [f"{source}:{cid}:gsas:v1", f"{source}:{cid}",
                            f"cod:{cid}", f"mp:{cid}:gsas:v1", f"mp:{cid}"]:
                text = _cache.get(key_fmt)
                if text:
                    cal_phase['cif_text'] = text
                    break

            # Override cell to NIST SRM 640g certified value for Si.
            # Only apply for SG 227 (Fd-3m / Si / Ge / diamond).
            _cal_sg = int(cal_phase.get('spacegroup_number', 0) or 0)
            if _cal_sg == 227:
                cal_phase['a'] = 5.431109
                cal_phase['b'] = 5.431109
                cal_phase['c'] = 5.431109
                cal_phase['alpha'] = 90.0
                cal_phase['beta']  = 90.0
                cal_phase['gamma'] = 90.0

            # Instrument from form (or auto-detect)
            _cal_instrument = form.get('instrument', '').strip().lower()
            if _cal_instrument in ('', 'auto'):
                _cal_instrument = None  # let run_calibration auto-detect

            cal_result = run_calibration(
                tt=data['tt'], y_obs=data['intensity'],
                sigma=data['sigma'],
                phase=cal_phase,
                wavelength=wavelength,
                tt_min=tt_min, tt_max=tt_max,
                n_bg_coeffs=int(form.get('n_bg_coeffs', 6)
                                if form.get('n_bg_coeffs', 'auto') != 'auto'
                                else 6),
                instrument=_cal_instrument,
                keep_workdir=True,
            )

            # Compute Si tick positions for the plot
            _cal_ticks = []
            try:
                from modules.xrd.crystallography import generate_reflections
                _cal_ticks_raw = generate_reflections(
                    cal_phase['a'], cal_phase['b'], cal_phase['c'],
                    cal_phase.get('alpha', 90), cal_phase.get('beta', 90),
                    cal_phase.get('gamma', 90),
                    'cubic', cal_phase.get('spacegroup_number', 227),
                    wavelength, tt_min, tt_max, hkl_max=12)
                _cal_ticks = [r[0] for r in _cal_ticks_raw]
            except Exception as _te:
                print(f"  Warning: tick positions: {_te}", flush=True)

            # Phase-only contribution = calc - background
            _y_calc = np.array(cal_result['y_calc'])
            _y_bg = np.array(cal_result['y_background'])
            _phase_only = np.maximum(_y_calc - _y_bg, 0).tolist()

            # Make a plot from the calibration result
            plot_result = {
                'tt': cal_result['tt'],
                'y_obs': cal_result['y_obs'],
                'y_calc': cal_result['y_calc'],
                'y_background': cal_result['y_background'],
                'residuals': cal_result['residuals'],
                'statistics': cal_result['statistics'],
                'phase_results': [{
                    'name': cal_phase.get('name', 'Standard'),
                    'spacegroup': cal_phase.get('spacegroup', 'Fd-3m'),
                    'weight_fraction_%': 100.0,
                    'tick_positions': _cal_ticks,
                }],
                'phase_patterns': [_phase_only],
                'wavelength': wavelength,
            }
            plot_path = os.path.join(out_dir, 'xrd_calibration.png')
            make_xrd_plot(plot_result, {
                'sample_id': sample_id,
                'wavelength_label': wl_label,
                'method': 'GSAS-II Calibration',
            }, plot_path)

            with open(plot_path, 'rb') as img:
                plot_b64 = base64.b64encode(img.read()).decode()

            return jsonify({
                'plot_b64':      plot_b64,
                'statistics':    cal_result['statistics'],
                'phase_results': [{
                    'name': cal_phase.get('name', 'Standard'),
                    'spacegroup': cal_phase.get('spacegroup', ''),
                    'weight_fraction_%': 100.0,
                    'instprm_params': cal_result['params'],
                    'lorentzian_term': cal_result.get('lorentzian_term', '?'),
                }],
                'zero_shift':    cal_result['params'].get('Zero', 0),
                'pymatgen_used': False,
                'method':        'GSAS-II Calibration',
                'instprm_path':  cal_result['instprm_path'],
                'output_dir':    out_dir,
            })

        # ── Normal refinement mode ─────────────────────────────────────
        print(f"  XRD backend call: method={form.get('method')} "
              f"instrument={form.get('instrument', 'auto')} "
              f"n_phases={len(phases)}", flush=True)
        result = xrd_processor.run(
            filepath   = upload_path,
            output_dir = out_dir,
            metadata   = {'sample_id': sample_id, 'notes': notes},
            params     = {
                'phases':           phases,
                'wavelength':       wavelength,
                'wavelength_label': wl_label,
                'tt_min':           tt_min,
                'tt_max':           tt_max,
                'n_bg_coeffs':      form.get('n_bg_coeffs', 'auto'),
                'max_outer':        MAX_OUTER,
                'method':           form.get('method', 'lebail'),
                'instprm_file':     instprm_file_path,
                'instrument':       form.get('instrument', 'auto'),
                # Verification mode (GSAS-II only): skip cell/Uiso/size
                # stages and refine only bg + scales + displacement + Y.
                # Use for first-pass tests when peak positions or widths
                # may be wrong — surfaces those problems instead of
                # absorbing them into more flexible parameters.
                'verification_mode': form.get('verification_mode', '').lower() == 'true',
                # Companion to verification_mode: when both are True,
                # Stage 6 cell refinement is enabled (Stages 3/4/4b/5
                # still skipped).  Use after a verification run shows
                # peaks land in the right place and you want the cell
                # to relax onto the experimental positions.
                'verify_refine_cell': form.get('verify_refine_cell', '').lower() == 'true',
                # GSAS-II phase isolation: when True, per-phase curves
                # come from GSAS-II's actual ycalc (each phase isolated
                # by zeroing other scales).  When False, per-phase
                # curves are reconstructed manually from the Python
                # reflection list (display-only, can disagree with
                # GSAS-II's internal Fc²).  Diagnostic toggle.
                'phase_isolation': form.get('phase_isolation', '').lower() == 'true',
                # Refine March-Dollase preferred orientation for
                # hexagonal/trigonal phases (e.g. WC along [001]) in
                # verification mode.  Forces preferred_orientation='auto'
                # internally so hex phases get MD enabled, then lets
                # Stage 6 refine the MD ratio alongside cell.
                'verify_refine_po':
                    form.get('verify_refine_po', '').lower() == 'true',
                # Swap position handle: refine Zero, fix DisplaceX/Y at 0.
                # Use when DisplaceY refuses to move from 0 because the
                # offset is actually a Zero miscalibration (the measured
                # instprm's Zero may not transfer cleanly to a different
                # sample mounting).  Overrides the measured-instprm rule
                # that locks Zero.
                'verify_use_zero_not_displace':
                    form.get('verify_use_zero_not_displace', '').lower() == 'true',
                # Branch B: post-Stage-6 enforce uniform cell scaling on
                # W2C.  After cell refines (a, b, c) freely, we override
                # them with (s·a₀, s·b₀, s·c₀) where s preserves the
                # refined volume.  Forces axis ratios back to mp-2034.
                # Diagnostic for whether W2C anisotropic refinement was
                # capturing real strain or just data-resolution noise.
                'verify_cell_uniform_w2c':
                    form.get('verify_cell_uniform_w2c', '').lower() == 'true',
                # Diagnostic: free X (Lorentzian strain).  Added to the
                # refine list in Stages 2 and 6.  Use to test whether
                # residual peak shape needs strain-like broadening on
                # top of the size broadening Y already captures.
                'verify_refine_x':
                    form.get('verify_refine_x', '').lower() == 'true',
                # Y controls (orthogonal to X / Zero handles).
                # fix_y dominates over Y-refinement when both set.
                # nonneg clamps post-Stage-6 Y < 0 to 0.
                'verify_fix_y':
                    form.get('verify_fix_y', '').lower() == 'true',
                'verify_y_fixed_value': (
                    float(form.get('verify_y_fixed_value', '0') or '0')
                    if form.get('verify_y_fixed_value', '').strip() else None),
                'verify_y_nonnegative':
                    form.get('verify_y_nonnegative', '').lower() == 'true',
                # Explicit structural toggles — override the verification_mode
                # default of "off" without leaving verification_mode.
                'verify_refine_uiso':
                    form.get('verify_refine_uiso', '').lower() == 'true',
                'verify_refine_size':
                    form.get('verify_refine_size', '').lower() == 'true',
                # Refine atom positions (XYZ).  Existing flag re-exposed.
                'refine_xyz':
                    form.get('refine_xyz', '').lower() == 'true',
                # Tick source: True = GSAS-II RefList, False = Python refs.
                'use_gsas_ref_ticks':
                    form.get('use_gsas_ref_ticks', '').lower() == 'true',
                # Fix WC PO: hold March-Dollase ratio at a user-specified
                # value during refinement (used to break MD ↔ Uiso
                # correlation in runs that combine PO and Uiso).
                'verify_fix_po':
                    form.get('verify_fix_po', '').lower() == 'true',
                'verify_po_fixed_value': (
                    float(form.get('verify_po_fixed_value', '0.905') or '0.905')
                    if form.get('verify_po_fixed_value', '').strip() else None),
                # Legacy per-phase HAP toggles (kept for back-compat).
                # Prefer phase_options (below) for new clients.
                'verify_refine_wc_size':
                    form.get('verify_refine_wc_size', '').lower() == 'true',
                'verify_refine_w2c_size':
                    form.get('verify_refine_w2c_size', '').lower() == 'true',
                'verify_refine_wc_mustrain':
                    form.get('verify_refine_wc_mustrain', '').lower() == 'true',
                'verify_refine_w2c_mustrain':
                    form.get('verify_refine_w2c_mustrain', '').lower() == 'true',
                # Generic per-phase refinement options.  JSON-serialized
                # list of dicts (one per phase, by index) with keys
                # refine_size, refine_mustrain, po_mode, po_value, po_axis.
                # Frontend builds this from the per-phase control cards.
                'phase_options':
                    (lambda _raw: (
                        __import__('json').loads(_raw)
                        if _raw and _raw.strip().startswith(('[', '{'))
                        else None
                    ))(form.get('phase_options', '')),
            }
        )
        print("=== xrd_processor.run returned ===", flush=True)

        with open(result['plot_path'], 'rb') as img:
            plot_b64 = base64.b64encode(img.read()).decode()

        print("=== /api/process_xrd DONE ===", flush=True)
        return jsonify({
            'plot_b64':      plot_b64,
            'statistics':    result['statistics'],
            'phase_results': result['phase_results'],
            'zero_shift':    result['zero_shift'],
            'displacement_um':    result.get('displacement_um'),
            'displacement_param': result.get('displacement_param'),
            'fit_warnings':  result.get('warnings', []),
            'pymatgen_used': result.get('pymatgen_used', False),
            'method':        result.get('method', 'Le Bail'),
            'summary_path':  result['summary_path'],
            'output_dir':    out_dir,
            # Raw arrays for interactive plotting — nested in result['result']
            'plot_data': {
                'tt':              result.get('result', {}).get('tt', []),
                'y_obs':           result.get('result', {}).get('y_obs', []),
                'y_calc':          result.get('result', {}).get('y_calc', []),
                'y_background':    result.get('result', {}).get('y_background', []),
                'residuals':       result.get('result', {}).get('residuals', []),
                'phase_patterns':  result.get('result', {}).get('phase_patterns', []),
            },
        })
    except (SystemExit, KeyboardInterrupt):
        # GSAS-II sometimes calls sys.exit() on fatal errors — catch it
        # so the Flask server stays alive.
        import traceback
        tb = traceback.format_exc()
        print(f"\n  !! XRD refinement crashed (SystemExit):\n{tb}", flush=True)
        return jsonify({
            'error': 'The refinement engine crashed unexpectedly. '
                     'Check the terminal for details. '
                     'If using GSAS-II, it may not be fully installed — '
                     'try Le Bail or Rietveld instead.',
            'trace': tb,
        }), 500
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"\n  !! XRD refinement error:\n{tb}", flush=True)
        return jsonify({'error': str(e), 'trace': tb}), 500

# ── Launch ────────────────────────────────────────────────────────────────────

def startup_url():
    target = (os.environ.get('CATALYSIS_TOOLKIT_START_PATH') or '/xrd').strip()
    if target.startswith(('http://', 'https://')):
        return target
    if not target.startswith('/'):
        target = '/' + target
    return f'http://localhost:5000{target}'


def open_browser():
    webbrowser.open(startup_url())

@app.route('/api/xrd/preview', methods=['POST'])
def xrd_preview():
    """
    Parse an uploaded XRD file and return data for live preview plot.
    Downsampled to ≤2000 points for fast rendering.
    """
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'No file'}), 400
        f = request.files['file']
        safe = re.sub(r'[^\w\-.]', '_', f.filename)
        path = os.path.join(UPLOAD_DIR, safe)
        f.save(path)

        data = xrd_processor.parse_xrd_file(path)
        tt   = [float(v) for v in data['tt']]
        y    = [float(v) for v in data['intensity']]

        MAX_PTS = 2000
        if len(tt) > MAX_PTS:
            step = len(tt) // MAX_PTS
            tt = tt[::step]
            y  = y[::step]

        return jsonify({
            'tt':        tt,
            'intensity': y,
            'tt_min':    round(min(tt), 2),
            'tt_max':    round(max(tt), 2),
            'y_max':     round(max(y), 1),
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    # Check GSAS-II status for startup banner
    try:
        from modules.xrd.gsasii_backend import is_available as _gsas_avail
        _gsas_status = 'ready' if _gsas_avail() else 'not installed'
    except Exception:
        _gsas_status = 'not installed'

    print("\n" + "━"*50)
    print("  Catalysis Data Toolkit")
    print(f"  pymatgen:          {'ready' if _pymatgen_ready else 'not installed'}")
    print(f"  GSAS-II:           {_gsas_status}")
    print(f"  Materials Project: {'configured' if MP_API_KEY else 'no API key'}")
    print(f"  CIF cache:         {_cache.stats()['entries']} entries "
          f"({_cache.stats()['size_mb']} MB)")
    print(f"  {startup_url()}")
    print("━"*50 + "\n")

    # Only open browser once — not in the reloader child process
    import os
    if os.environ.get("WERKZEUG_RUN_MAIN") != "true":
        Timer(1.2, open_browser).start()

    app.run(debug=True, port=5000, use_reloader=False, threaded=True)

