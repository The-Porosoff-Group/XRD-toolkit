"""
modules/xrd/mp_api.py
Materials Project database interface — new API (api.materialsproject.org).

Your key (from next-gen.materialsproject.org/api) works ONLY with:
  https://api.materialsproject.org
  Header: X-API-KEY: <key>

The v2 legacy API (materialsproject.org/rest/v2) requires a SEPARATE legacy
key from legacy.materialsproject.org/open — do not mix them up.

The current summary endpoint accepts ``structure`` in ``_fields``.  Structure
records are normalized to conventional, symmetry-bearing CIF text before they
enter the rest of the toolkit.
"""

import math, os, re, requests, json
from .crystallography import parse_cif, conventionalize_phase_cell
from .cod_api import infer_system, _sf
from .cif_cache import normalize_mp_id

MP_SUMMARY = "https://api.materialsproject.org/materials/summary/"
TIMEOUT    = 15

_VALID_ELEMENTS = {
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


def _normalize_formula_case(formula):
    formula = (formula or '').strip().replace(' ', '')
    if not formula or re.search(r'[A-Z]', formula):
        return formula
    out = []
    i = 0
    while i < len(formula):
        if formula[i].isdigit():
            out.append(formula[i])
            i += 1
            continue
        two = formula[i:i+2].capitalize()
        one = formula[i].upper()
        if i + 1 < len(formula) and two in _VALID_ELEMENTS:
            out.append(two)
            i += 2
        elif one in _VALID_ELEMENTS:
            out.append(one)
            i += 1
        else:
            out.append(formula[i].upper())
            i += 1
    return ''.join(out)

# ─────────────────────────────────────────────────────────────────────────────
# Local CIF fixtures override
# ─────────────────────────────────────────────────────────────────────────────
# Some Materials Project entries import incorrectly into GSAS-II when the raw
# structure JSON is round-tripped through pymatgen's CifWriter (e.g. mp-2034
# W2C used to land as P1/full-cell, blowing up the cell DoF count).  Audited
# canonical CIFs live in fixtures/ at the toolkit root.  When fetch_cif() sees
# one of these mp_ids, it substitutes the fixture text for the round-tripped
# CIF, but keeps everything else (formula, symmetry metadata, etc.) from MP.
#
# To add a new fixture: drop the .cif into fixtures/ and add an entry below.
_FIXTURE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    'fixtures')
_LOCAL_FIXTURES = {
    'mp-2034': 'w2c_pbcn_mp_2034.cif',   # W2C Pbcn — see CIF-Audit_v1.md
    'mp-129': 'mo_metal_bcc_im3m.cif',   # Mo bcc Im-3m — use RT conventional cell
    'mp-1894': 'wc_p-6m2_mp_1894.cif',
    'mp-33065': 'w2c_pnnm_mp_33065.cif',
    'mp-684989': 'w9c4_r32_mp_684989.cif',
    'mp-567397': 'w2c_p-31m_mp_567397.cif',
    'mp-1008625': 'w2c_p-3m1_mp_1008625.cif',
    'mp-13136': 'gamma_wc1x_fm3m.cif',
    'mp-1552': 'mo2c_pbcn_mp_1552.cif',
    'mp-2305': 'moc_p-6m2_mp_2305.cif',
    'mp-1221498': 'mo2c_p-3m1_mp_1221498.cif',
    'mp-1221473': 'mo3c2_p-3m1_mp_1221473.cif',
    'mp-2746': 'gamma_moc1x_fm3m.cif',
}

_LOCAL_FIXTURE_METADATA = {
    'mp-2034': {'formula': 'W2C', 'spacegroup': 'Pbcn',
                'spacegroup_number': 60, 'system': 'orthorhombic', 'Z': 4},
    'mp-129': {'formula': 'Mo', 'spacegroup': 'Im-3m',
               'spacegroup_number': 229, 'system': 'cubic', 'Z': 2},
    'mp-1894': {'formula': 'WC', 'spacegroup': 'P-6m2',
                'spacegroup_number': 187, 'system': 'hexagonal', 'Z': 1},
    'mp-33065': {'formula': 'W2C', 'spacegroup': 'Pnnm',
                 'spacegroup_number': 58, 'system': 'orthorhombic', 'Z': 2},
    'mp-684989': {'formula': 'W9C4', 'spacegroup': 'R32',
                  'spacegroup_number': 155, 'system': 'trigonal', 'Z': 6},
    'mp-567397': {'formula': 'W2C', 'spacegroup': 'P-31m',
                  'spacegroup_number': 162, 'system': 'trigonal', 'Z': 3},
    'mp-1008625': {'formula': 'W2C', 'spacegroup': 'P-3m1',
                   'spacegroup_number': 164, 'system': 'trigonal', 'Z': 1},
    'mp-13136': {'formula': 'WC1-x', 'spacegroup': 'Fm-3m',
                 'spacegroup_number': 225, 'system': 'cubic', 'Z': 4,
                 'gamma_wc1x': True},
    'mp-1552': {'formula': 'Mo2C', 'spacegroup': 'Pbcn',
                'spacegroup_number': 60, 'system': 'orthorhombic', 'Z': 4},
    'mp-2305': {'formula': 'MoC', 'spacegroup': 'P-6m2',
                'spacegroup_number': 187, 'system': 'hexagonal', 'Z': 1},
    'mp-1221498': {'formula': 'Mo2C', 'spacegroup': 'P-3m1',
                   'spacegroup_number': 164, 'system': 'trigonal', 'Z': 1},
    'mp-1221473': {'formula': 'Mo3C2', 'spacegroup': 'P-3m1',
                   'spacegroup_number': 164, 'system': 'trigonal', 'Z': 1},
    'mp-2746': {'formula': 'MoC1-x', 'spacegroup': 'Fm-3m',
                'spacegroup_number': 225, 'system': 'cubic', 'Z': 4,
                'gamma_wc1x': True},
}


def _fixture_cif_for(mp_id, purpose=None):
    """Return fixture CIF text, optionally enforcing its declared role.

    ``purpose='normal_import'`` permits only audited fixtures whose manifest
    includes ``normal_import`` and does not mark ``normal_import_safe`` false.
    Raw/P1 fixtures remain accessible to audit and regression callers that do
    not request a production purpose.
    """
    mp_id = normalize_mp_id(mp_id)
    fname = _LOCAL_FIXTURES.get(mp_id)
    if not fname:
        return None
    if purpose == 'normal_import':
        record = _fixture_record_for(mp_id)
        intended = set(record.get('intended_use') or [])
        if ('normal_import' not in intended
                or record.get('normal_import_safe', True) is False):
            return None
    path = os.path.join(_FIXTURE_DIR, fname)
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as f:
            return f.read()
    except Exception:
        return None


def _fixture_metadata_for(mp_id):
    """Return curated metadata for a local fixture, if any."""
    return dict(_LOCAL_FIXTURE_METADATA.get(normalize_mp_id(mp_id), {}))


def _fixture_manifest():
    """Load fixture provenance/role metadata, if present."""
    path = os.path.join(_FIXTURE_DIR, 'fixture_manifest.json')
    try:
        with open(path, encoding='utf-8') as f:
            return json.load(f).get('fixtures', {})
    except Exception:
        return {}


def _fixture_record_for(mp_id):
    """Return manifest record for the fixture mapped to an mp-id."""
    fname = _LOCAL_FIXTURES.get(normalize_mp_id(mp_id))
    if not fname:
        return {}
    return dict(_fixture_manifest().get(fname, {}))


def _apply_fixture_record(result, mp_id):
    """Attach fixture role metadata and warnings to a phase dict."""
    rec = _fixture_record_for(mp_id)
    if not rec:
        return result
    result['fixture_cell_setting'] = rec.get('cell_setting')
    result['fixture_intended_use'] = rec.get('intended_use') or []
    result['fixture_normal_import_safe'] = rec.get('normal_import_safe', True)
    if ('normal_import' in result['fixture_intended_use']
            and result['fixture_normal_import_safe']):
        result['cif_preparation_policy'] = 'audited_normal_fixture'
    if rec.get('notes'):
        result['fixture_notes'] = rec.get('notes')
    if rec.get('normal_import_safe') is False:
        result['fixture_warning'] = (
            'This local fixture is a raw/source CIF, not a validated '
            'canonical conventional CIF. Preview may be useful, but GSAS-II '
            'refinement should be checked carefully.'
        )
    return result


# ─────────────────────────────────────────────────────────────────────────────
# SEARCH
# ─────────────────────────────────────────────────────────────────────────────

def _get(params, api_key):
    """Core GET to the summary endpoint."""
    headers = {"X-API-KEY": api_key, "Accept": "application/json"}
    # '_fields' is the correct parameter for field selection in the new MP REST API
    # 'fields' (without underscore) causes a 400; default returns only material_id
    params["_fields"]    = ("material_id,formula_pretty,symmetry,"
                             "energy_above_hull,theoretical,nsites,nelements,"
                             "volume,density,structure")
    params["deprecated"] = "false"
    params["_limit"]     = params.get("_limit", 50)

    resp = requests.get(MP_SUMMARY, headers=headers,
                         params=params, timeout=TIMEOUT)

    if resp.status_code == 403:
        return {"error": "Materials Project API key invalid or expired. "
                         "Check config.yaml — make sure it is the key from "
                         "next-gen.materialsproject.org/api, not legacy."}
    if resp.status_code == 400:
        return {"error": f"Bad request to Materials Project API. "
                         f"Details: {resp.text[:200]}"}
    resp.raise_for_status()
    return resp.json().get("data", [])


def search_by_elements(elements, api_key, strict=True,
                        max_results=50, sort_by="formula"):
    if not api_key:
        return {"error": "No Materials Project API key. Add to config.yaml."}
    elements = [e.strip().capitalize() for e in elements if e.strip()]
    if not elements:
        return {"error": "No elements provided."}
    try:
        if strict:
            params = {"chemsys": "-".join(sorted(elements)), "_limit": max_results}
        else:
            params = {"elements": ",".join(elements), "_limit": max_results}
        data = _get(params, api_key)
        if isinstance(data, dict) and "error" in data:
            return data
        return _sort(_parse(data), sort_by)
    except requests.exceptions.ConnectionError:
        return {"error": "Cannot reach Materials Project. Check internet."}
    except requests.exceptions.Timeout:
        return {"error": "Materials Project search timed out."}
    except Exception as e:
        return {"error": f"Materials Project search error: {e}"}


def search_by_formula(formula, api_key, max_results=50, sort_by="formula"):
    if not api_key:
        return {"error": "No Materials Project API key. Add to config.yaml."}
    formula = _normalize_formula_case(formula)
    if not formula:
        return []
    try:
        # Try exact formula first
        data = _get({"formula": formula, "_limit": max_results}, api_key)
        if isinstance(data, dict) and "error" in data:
            return data
        results = _parse(data)
        if results:
            return _sort(results, sort_by)

        # Fallback: chemsys from elements in formula
        elements = list(dict.fromkeys(re.findall(r"[A-Z][a-z]?", formula)))
        if elements:
            data2 = _get({"chemsys": "-".join(sorted(elements)),
                           "_limit": max_results}, api_key)
            if isinstance(data2, dict) and "error" in data2:
                return data2
            return _sort(_parse(data2), sort_by)
        return []
    except requests.exceptions.ConnectionError:
        return {"error": "Cannot reach Materials Project. Check internet."}
    except requests.exceptions.Timeout:
        return {"error": "Materials Project search timed out."}
    except Exception as e:
        return {"error": f"Materials Project search error: {e}"}


def search_by_name(name, api_key, max_results=50, sort_by="formula"):
    """
    MP has no free-text search. Routes by input type:
    - Looks like a formula (W2C, WC, Mo2C) → formula search
    - Contains known element symbols (W, Mo, Fe) → chemsys search
    - Plain English names → try to map common words to elements, then chemsys
    """
    if not api_key:
        return {"error": "No Materials Project API key. Add to config.yaml."}
    name = name.strip()
    if not name:
        return []

    # Common element name → symbol mappings
    _NAME_MAP = {
        "tungsten": "W", "molybdenum": "Mo", "iron": "Fe", "carbon": "C",
        "nitrogen": "N", "oxygen": "O", "silicon": "Si", "nickel": "Ni",
        "cobalt": "Co", "copper": "Cu", "chromium": "Cr", "vanadium": "V",
        "titanium": "Ti", "zirconium": "Zr", "hafnium": "Hf", "niobium": "Nb",
        "tantalum": "Ta", "rhenium": "Re", "ruthenium": "Ru", "palladium": "Pd",
        "platinum": "Pt", "gold": "Au", "silver": "Ag", "aluminium": "Al",
        "aluminum": "Al", "manganese": "Mn", "zinc": "Zn", "tin": "Sn",
        "lead": "Pb", "sulfur": "S", "phosphorus": "P", "boron": "B",
        "carbide": "C", "nitride": "N", "oxide": "O", "silicide": "Si",
    }

    # If it looks like a formula (starts uppercase, only letters/digits, no spaces)
    _compact_name = name.replace(" ", "")
    _maybe_formula = _normalize_formula_case(_compact_name)
    if (len(_compact_name) <= 8
            and re.match(r"^[A-Z][a-zA-Z0-9]*$", _maybe_formula)):
        return search_by_formula(_maybe_formula, api_key, max_results, sort_by)

    # Try to extract element symbols — first from capitalised tokens (e.g. "W C Mo")
    words = name.replace("-", " ").split()
    elements = []
    for word in words:
        w = word.strip("(),.")
        # Direct element symbol match (1-2 chars, starts uppercase)
        if re.match(r"^[A-Z][a-z]?$", w):
            elements.append(w)
        # English name lookup
        elif w.lower() in _NAME_MAP:
            el = _NAME_MAP[w.lower()]
            if el not in elements:
                elements.append(el)

    if elements:
        return search_by_elements(elements, api_key, strict=True,
                                   max_results=max_results, sort_by=sort_by)
    return []


# ─────────────────────────────────────────────────────────────────────────────
# PARSE
# ─────────────────────────────────────────────────────────────────────────────

def _parse(entries):
    results = []
    for e in entries:
        try:
            api_mp_id = str(e.get("material_id", ""))
            mp_id = normalize_mp_id(api_mp_id)
            formula = str(e.get("formula_pretty", ""))
            sym     = e.get("symmetry") or {}
            sg_sym  = str(sym.get("symbol", ""))
            sg_num  = int(sym.get("number") or 1)
            cs      = (sym.get("crystal_system") or "").lower()
            system  = cs or infer_system(sg_num)
            e_hull  = float(e.get("energy_above_hull") or 0)

            # Extract cell params from structure.lattice if present
            a = b = c = al = be = ga = None
            struct = e.get("structure") or {}
            lattice = struct.get("lattice") or {}
            if lattice:
                a  = _sf(lattice.get("a"))
                b  = _sf(lattice.get("b"))
                c  = _sf(lattice.get("c"))
                al = _sf(lattice.get("alpha"), 90.0)
                be = _sf(lattice.get("beta"),  90.0)
                ga = _sf(lattice.get("gamma"), 90.0)

            result = conventionalize_phase_cell({
                "mp_id":             mp_id,
                "mp_api_id":         api_mp_id,
                "cod_id":            mp_id,
                "formula":           formula,
                "name":              formula,
                "spacegroup":        sg_sym,
                "spacegroup_number": sg_num,
                "system":            system or "triclinic",
                "a": a, "b": b, "c": c,
                "alpha": al or 90.0,
                "beta":  be or 90.0,
                "gamma": ga or 90.0,
                "stability":         _stab(e_hull),
                "e_above_hull":      round(e_hull, 4),
                "theoretical":       bool(e.get("theoretical", True)),
                "year":              "DFT",
                "authors":           "Materials Project",
                "journal":           "Comp.",
                "source":            "mp",
            })

            fixture_text = _fixture_cif_for(
                mp_id, purpose='normal_import')
            if fixture_text:
                fixture = parse_cif(fixture_text)
                fixture_meta = _fixture_metadata_for(mp_id)
                fixture_sg = int(fixture.get('spacegroup_number') or 1)
                if fixture_sg > 1:
                    for key in ('formula', 'a', 'b', 'c', 'alpha', 'beta',
                                'gamma', 'spacegroup_number', 'system', 'Z'):
                        if fixture.get(key) not in (None, ''):
                            result[key] = fixture[key]
                    result['spacegroup'] = (
                        fixture.get('spacegroup')
                        or fixture.get('spacegroup_name')
                        or result.get('spacegroup')
                    )
                    result['name'] = fixture.get('formula') or result.get('name')
                for key, value in fixture_meta.items():
                    if value not in (None, ''):
                        result[key] = value
                result = conventionalize_phase_cell(result)
                if fixture_meta.get('gamma_wc1x'):
                    result['name'] = fixture_meta.get('formula') or result.get('name')
                result['_cif_text'] = fixture_text
                result = _apply_fixture_record(result, mp_id)

                try:
                    for el, x, y, z, occ in fixture.get('sites') or []:
                        if (str(el).upper() == 'C'
                                and abs(float(x) % 1.0 - 0.5) < 1e-4
                                and abs(float(y) % 1.0 - 0.5) < 1e-4
                                and abs(float(z) % 1.0 - 0.5) < 1e-4):
                            result['gamma_c_occupancy'] = float(occ)
                            result['gamma_vacancy_x'] = max(
                                0.0, min(1.0, 1.0 - float(occ)))
                            result['gamma_wc1x'] = True
                            break
                except Exception:
                    pass

            results.append(result)
        except Exception:
            continue
    return results


def _stab(e):
    if e < 0.001:  return "stable (on hull)"
    elif e < 0.05: return f"near-stable (+{e*1000:.0f} meV/atom)"
    elif e < 0.15: return f"metastable (+{e*1000:.0f} meV/atom)"
    else:          return f"unstable (+{e*1000:.0f} meV/atom)"


def _sort(results, sort_by):
    if sort_by == "year_desc":
        results.sort(key=lambda r: r.get("e_above_hull", 99))
    elif sort_by == "cell_a":
        results.sort(key=lambda r: r.get("a") or 999)
    elif sort_by == "spacegroup":
        results.sort(key=lambda r: r.get("spacegroup_number", 999))
    else:
        results.sort(key=lambda r: r.get("formula", ""))
    return results


# ─────────────────────────────────────────────────────────────────────────────
# CIF DOWNLOAD
# ─────────────────────────────────────────────────────────────────────────────

def fetch_cif(mp_id, api_key):
    """
    Fetch structure for a Materials Project entry and convert to CIF.
    The new API has no dedicated /cif endpoint — we request 'structure'
    from the summary endpoint, parse it with pymatgen, and use spglib to write
    a conventional symmetry-bearing CIF.
    """
    requested_mp_id = str(mp_id)
    mp_id = normalize_mp_id(requested_mp_id)
    fixture_text = _fixture_cif_for(mp_id, purpose='normal_import')
    if fixture_text:
        print(f"  fetch_cif: using local fixture for {mp_id} "
              f"(no MP API call needed)", flush=True)
        parsed = parse_cif(fixture_text)
        parsed.update(_fixture_metadata_for(mp_id))
        parsed = conventionalize_phase_cell(parsed)
        parsed.update({"mp_id": mp_id, "cod_id": mp_id,
                       "cif_text": fixture_text, "source": "mp"})
        parsed = _apply_fixture_record(parsed, mp_id)
        return parsed

    if not api_key:
        raise ValueError("No Materials Project API key configured.")

    headers = {"X-API-KEY": api_key, "Accept": "application/json"}

    # Request structure from summary endpoint
    resp = requests.get(
        MP_SUMMARY,
        headers=headers,
        params={
            "material_ids": mp_id,
            "_fields":      "material_id,formula_pretty,symmetry,structure",
            "deprecated":   "false",
            "_limit":       1,
        },
        timeout=TIMEOUT,
    )
    if resp.status_code == 403:
        raise PermissionError("API key invalid or expired.")
    resp.raise_for_status()

    data    = resp.json().get("data", [])
    if not data:
        raise RuntimeError(f"No entry found for {mp_id}")

    entry   = data[0]
    api_mp_id = str(entry.get('material_id') or requested_mp_id)
    struct  = entry.get("structure")
    formula = entry.get("formula_pretty", "")
    sym     = entry.get("symmetry") or {}

    if not struct:
        raise RuntimeError(f"No structure data returned for {mp_id} "
                           f"— 'structure' may not be available via raw REST")

    # Convert the pymatgen structure dict to conventional CIF text.
    cif_text = _structure_dict_to_cif(struct, mp_id, formula, sym)

    parsed = parse_cif(cif_text)
    _source_cif_sg = int(parsed.get('spacegroup_number') or 1)
    parsed.update({"mp_id": mp_id, "mp_api_id": api_mp_id, "cod_id": mp_id,
                   "formula": formula, "cif_text": cif_text, "source": "mp"})
    parsed['source_cif_spacegroup_number'] = _source_cif_sg
    parsed['cif_preparation_policy'] = (
        'mp_conventional_cif'
        if _source_cif_sg > 1
        else 'mp_p1_full_cell_fallback')

    # Merge the authoritative MP symmetry metadata.  The direct-spglib path
    # writes matching symmetry into the CIF; the P1 emergency fallback does
    # not, so this also keeps result metadata informative in that case.
    if sym:
        if sym.get('number'):
            parsed['spacegroup_number'] = int(sym['number'])
        if sym.get('symbol'):
            parsed['spacegroup'] = sym['symbol']
            parsed['spacegroup_name'] = sym['symbol']
        if sym.get('crystal_system'):
            parsed['system'] = sym['crystal_system'].lower()

    return parsed


def _element_symbol(specie):
    symbol = getattr(specie, 'symbol', None)
    if symbol:
        return str(symbol)
    element = getattr(specie, 'element', None)
    if element is not None and getattr(element, 'symbol', None):
        return str(element.symbol)
    match = re.match(r'([A-Z][a-z]?)', str(specie))
    return match.group(1) if match else 'X'


def _symop_component(coeffs, translation):
    """Format one spglib affine component as a CIF xyz expression."""
    from fractions import Fraction
    terms = []
    for coefficient, variable in zip(coeffs, ('x', 'y', 'z')):
        coefficient = int(coefficient)
        if not coefficient:
            continue
        sign = '-' if coefficient < 0 else '+'
        magnitude = abs(coefficient)
        body = variable if magnitude == 1 else f'{magnitude}{variable}'
        terms.append((sign, body))
    fraction = Fraction(float(translation) % 1.0).limit_denominator(24)
    if fraction:
        terms.append(('+', (str(fraction.numerator) if fraction.denominator == 1
                            else f'{fraction.numerator}/{fraction.denominator}')))
    if not terms:
        return '0'
    first_sign, first_body = terms[0]
    text = ('-' if first_sign == '-' else '') + first_body
    for sign, body in terms[1:]:
        text += sign + body
    return text


def _spglib_conventional_cif(struct_dict, mp_id, formula, sym):
    """Build a matching conventional/asymmetric CIF without pymatgen SGA.

    ``SpacegroupAnalyzer.get_conventional_standard_structure`` and
    ``CifWriter(symprec=...)`` can terminate the interpreter for valid oxide
    structures in the bundled Windows environment.  Direct spglib calls are
    stable.  They also provide the conventional cell, equivalence classes,
    and exact symmetry operations in one internally consistent setting.
    """
    import numpy as np
    import spglib
    from pymatgen.core import Structure

    structure = Structure.from_dict(struct_dict)
    signatures = []
    signature_ids = {}
    site_types = []
    for site in structure:
        signature = tuple(sorted(
            (_element_symbol(specie), round(float(occupancy), 8))
            for specie, occupancy in site.species.items()))
        if signature not in signature_ids:
            signature_ids[signature] = len(signatures) + 1
            signatures.append(signature)
        site_types.append(signature_ids[signature])

    cell = (np.asarray(structure.lattice.matrix, dtype=float),
            np.asarray(structure.frac_coords, dtype=float), site_types)
    expected_sg = int((sym or {}).get('number') or 1)

    for symprec in (0.01, 0.05, 0.1, 0.2):
        standardized = spglib.standardize_cell(
            cell, to_primitive=False, no_idealize=False, symprec=symprec)
        if standardized is None:
            continue
        dataset = spglib.get_symmetry_dataset(standardized, symprec=symprec)
        if dataset is None:
            continue
        detected_sg = int(dataset.number)
        if expected_sg > 1 and detected_sg != expected_sg:
            continue

        lattice, positions, standardized_types = standardized
        lengths = [float(np.linalg.norm(vector)) for vector in lattice]

        def angle(first, second):
            denominator = np.linalg.norm(first) * np.linalg.norm(second)
            cosine = float(np.dot(first, second) / denominator)
            return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))

        alpha = angle(lattice[1], lattice[2])
        beta = angle(lattice[0], lattice[2])
        gamma = angle(lattice[0], lattice[1])
        volume = abs(float(
            lattice[0][0] * (lattice[1][1] * lattice[2][2]
                             - lattice[1][2] * lattice[2][1])
            - lattice[0][1] * (lattice[1][0] * lattice[2][2]
                               - lattice[1][2] * lattice[2][0])
            + lattice[0][2] * (lattice[1][0] * lattice[2][1]
                               - lattice[1][1] * lattice[2][0])))
        sg_symbol = str(dataset.international).strip() or str(
            (sym or {}).get('symbol') or 'P 1')

        representative_indices = []
        seen_orbits = set()
        for index, orbit in enumerate(dataset.equivalent_atoms):
            orbit = int(orbit)
            if orbit not in seen_orbits:
                seen_orbits.add(orbit)
                representative_indices.append(index)

        atom_lines = []
        label_counts = {}
        for index in representative_indices:
            signature = signatures[int(standardized_types[index]) - 1]
            x, y, z = (float(value) % 1.0 for value in positions[index])
            for element, occupancy in signature:
                label_counts[element] = label_counts.get(element, 0) + 1
                atom_lines.append(
                    f" {element}{label_counts[element]:<6} {element:<4} "
                    f"{x:.8f} {y:.8f} {z:.8f} {occupancy:.6f}")

        symop_lines = []
        seen_operations = set()
        for rotation, translation in zip(dataset.rotations,
                                         dataset.translations):
            operation = ','.join(
                _symop_component(rotation[row], translation[row])
                for row in range(3))
            if operation not in seen_operations:
                seen_operations.add(operation)
                symop_lines.append(f" '{operation}'")

        if not atom_lines or not symop_lines:
            continue
        print(f"  MP CIF ({mp_id}): spglib standardized to SG "
              f"{detected_sg}, {len(standardized_types)} full-cell sites, "
              f"{len(representative_indices)} asymmetric site(s) "
              f"(symprec={symprec}).", flush=True)
        symop_text = '\n'.join(symop_lines)
        atom_text = '\n'.join(atom_lines)
        return f"""data_{mp_id}
_audit_creation_method           'Catalysis Toolkit direct spglib standardization'
_cell_length_a                   {lengths[0]:.8f}
_cell_length_b                   {lengths[1]:.8f}
_cell_length_c                   {lengths[2]:.8f}
_cell_angle_alpha                {alpha:.6f}
_cell_angle_beta                 {beta:.6f}
_cell_angle_gamma                {gamma:.6f}
_cell_volume                     {volume:.8f}
_symmetry_space_group_name_H-M   '{sg_symbol}'
_symmetry_Int_Tables_number      {detected_sg}
_chemical_formula_sum            '{formula}'
loop_
 _space_group_symop_operation_xyz
{symop_text}
loop_
 _atom_site_label
 _atom_site_type_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
{atom_text}
"""
    return None


def _structure_dict_to_cif(struct_dict, mp_id, formula, sym):
    """Convert a Materials Project structure JSON dict to a safe CIF.

    The primary route standardizes the cell and reduces the asymmetric unit
    with direct spglib data, rejecting any result whose detected space group
    differs from the MP summary record.  A plain P1 CifWriter dump is the safe
    fallback; space-group headers are never patched onto full-cell sites.
    """
    try:
        standardized = _spglib_conventional_cif(
            struct_dict, mp_id, formula, sym)
        if standardized:
            return standardized
    except BaseException as exc:
        print(f"  MP CIF ({mp_id}): direct spglib standardization failed "
              f"({exc}); using a self-consistent P1 fallback.", flush=True)

    try:
        from pymatgen.core import Structure
        from pymatgen.io.cif import CifWriter
        structure = Structure.from_dict(struct_dict)
        cif_text = str(CifWriter(structure))
        print(f"  MP CIF ({mp_id}): wrote P1 fallback with "
              f"{len(structure)} source-cell sites.", flush=True)
        return cif_text
    except BaseException:
        pass

    # Fallback: build a minimal CIF from the lattice dict
    lattice = struct_dict.get("lattice", {})
    a  = lattice.get("a",  4.0)
    b  = lattice.get("b",  a)
    c  = lattice.get("c",  a)
    al = lattice.get("alpha",  90.0)
    be = lattice.get("beta",   90.0)
    ga = lattice.get("gamma",  90.0)
    # This fallback writes every source-cell site, so it must remain P1.
    # Declaring the MP summary SG here would make downstream readers expand
    # the already-complete site list a second time.
    sg_num = 1
    sg_sym = "P 1"

    sites = struct_dict.get("sites", [])
    atom_lines = ""
    for i, site in enumerate(sites):
        species = site.get("species", [{}])
        abc = site.get("abc", [0, 0, 0])
        for j, component in enumerate(species):
            sp = component.get("element", "X")
            occupancy = float(component.get("occu", 1.0) or 1.0)
            atom_lines += (f" {sp}{i+1}_{j+1:<4} {sp:<4} "
                           f"{abc[0]:.6f} {abc[1]:.6f} {abc[2]:.6f} "
                           f"{occupancy:.6f}\n")

    return f"""data_{mp_id}
_cell_length_a                  {a:.6f}
_cell_length_b                  {b:.6f}
_cell_length_c                  {c:.6f}
_cell_angle_alpha               {al:.4f}
_cell_angle_beta                {be:.4f}
_cell_angle_gamma               {ga:.4f}
_symmetry_space_group_name_H-M  '{sg_sym}'
_symmetry_Int_Tables_number     {sg_num}
_chemical_formula_sum           '{formula}'
loop_
 _atom_site_label
 _atom_site_type_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
{atom_lines}"""


# ─────────────────────────────────────────────────────────────────────────────
# KEY VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

def validate_api_key(api_key):
    if not api_key or len(api_key) < 10:
        return False, "API key too short or missing."
    try:
        data = _get({"formula": "W", "_limit": 1}, api_key)
        if isinstance(data, dict) and "error" in data:
            return False, data["error"]
        return True, "API key valid."
    except requests.exceptions.ConnectionError:
        return False, "Cannot reach Materials Project (no internet?)."
    except Exception as e:
        return False, f"Validation error: {e}"
