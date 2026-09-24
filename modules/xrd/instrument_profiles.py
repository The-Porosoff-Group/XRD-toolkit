"""Built-in geometry presets and user-owned, local GSAS-II profiles.

User profiles live outside the source-controlled instrument presets. Each save
creates a new id, so neither bundled files nor earlier calibrations are replaced.
"""
import json
import math
import os
from pathlib import Path
import uuid
import warnings


TOOLKIT_ROOT = Path(__file__).resolve().parents[2]
LOCAL_INSTRUMENT_DIR = TOOLKIT_ROOT / 'local_instruments'
DEFAULT_INSTRUMENT = 'generic_flat_plate'


def geometry_profile(geometry):
    if geometry not in ('bragg_brentano', 'capillary'):
        raise ValueError('Choose flat plate / Bragg-Brentano or capillary / transmission.')
    flat = geometry == 'bragg_brentano'
    return {
        'geometry': geometry, 'displacement_param': 'Shift' if flat else 'DisplaceY',
        'zero_seed': 0.0, 'polariz': 0.5, 'sh_l': 0.002,
        'preferred_orientation_default': 'auto' if flat else 'off',
        'sigma_inflation_K': 1.0,
        'notes': 'Generic starting parameters; no measured instrument profile.',
    }


INSTRUMENT_PROFILES = {
    'generic_flat_plate': dict(geometry_profile('bragg_brentano'),
        label='Default / generic flat plate (Bragg-Brentano)'),
    'generic_capillary': dict(geometry_profile('capillary'),
        label='Generic capillary / transmission'),
    'benchtop_cu': dict(geometry_profile('bragg_brentano'),
        label='Benchtop Cu — flat plate (Si 640g)',
        zero_seed=-0.04099, sh_l=0.06399,
        instprm_filename='benchtop_Cu_Si640g.instprm',
        wavelength=1.540593, spectrum='cu_doublet',
        calibration_range=[20.0, 90.0],
        notes='Measured Si 640g calibration for the supplied benchtop instrument, '
              '2026-09-24; Rwp 8.50%. Make/model not recorded. Use only with '
              'the same instrument and optics; see docs/benchtop_calibration.md.'),
    'smartlab': dict(geometry_profile('bragg_brentano'),
        label='Rigaku SmartLab (BB)', zero_seed=-0.027, polariz=0.7,
        calibration_allow_x=False, calibration_allow_y=False,
        calibration_refine_sh_l=False, calibration_fixed_sh_l=0.002,
        calibration_u_min_rwp_gain=0.0, calibration_v_min_rwp_gain=0.0,
        instprm_filename='smartlab_Si640g.instprm',
        notes='Bundled SmartLab profile; use only with its matching configuration.'),
    'synergy_s': dict(geometry_profile('capillary'),
        label='Synergy-S (capillary)', zero_seed=-0.25, sigma_inflation_K=5.0,
        instprm_filename='synergy_s_Si640g.instprm',
        notes='Bundled Synergy-S profile; use only with its matching configuration.'),
}


def parse_instprm(content):
    """Validate a single-bank constant-wavelength X-ray file before using it."""
    if isinstance(content, bytes):
        try:
            content = content.decode('utf-8-sig')
        except UnicodeDecodeError as exc:
            raise ValueError('The .instprm file must be UTF-8 text.') from exc
    if not isinstance(content, str) or not content.strip() or len(content) > 100000:
        raise ValueError('Provide a nonempty .instprm file smaller than 100 KB.')
    values = {}
    for line in content.splitlines():
        if line.lstrip().startswith('#'):
            continue
        for item in line.split(';'):
            if ':' not in item:
                continue
            key, value = (part.strip() for part in item.split(':', 1))
            if key in values:
                raise ValueError(f'Duplicate {key}: choose a single-bank .instprm file.')
            values[key] = value
    if values.get('Type') != 'PXC':
        raise ValueError('Only GSAS-II constant-wavelength X-ray (Type:PXC) .instprm files are supported.')
    required = ['Zero', 'U', 'V', 'W', 'X', 'Y', 'SH/L']
    doublet = 'Lam1' in values or 'Lam2' in values
    required += ['Lam1', 'Lam2', 'I(L2)/I(L1)'] if doublet else ['Lam']
    if doublet and 'Lam' in values:
        raise ValueError('Use either Lam or Lam1/Lam2, not both.')
    for key in set(required + [key for key in ('Polariz.', 'Azimuth', 'Z') if key in values]):
        try:
            number = float(values[key])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f'Missing or invalid instrument parameter: {key}.') from exc
        if not math.isfinite(number):
            raise ValueError(f'Instrument parameter {key} must be finite.')
        values[key] = number
    for key in (['Lam1', 'Lam2'] if doublet else ['Lam']):
        if values[key] <= 0:
            raise ValueError('Instrument wavelengths must be positive.')
    for key in ('X', 'Y', 'SH/L', 'I(L2)/I(L1)'):
        if key in values and values[key] < 0:
            raise ValueError(f'Instrument parameter {key} cannot be negative.')
    if not 0 <= values.get('Polariz.', 0.5) <= 1:
        raise ValueError('Polarization must be between 0 and 1.')
    return values


def get_instrument_profiles():
    profiles = {key: dict(value) for key, value in INSTRUMENT_PROFILES.items()}
    for metadata in sorted(LOCAL_INSTRUMENT_DIR.glob('local_*.json')):
        try:
            saved = json.loads(metadata.read_text(encoding='utf-8'))
            key = metadata.stem
            profile = geometry_profile(saved['geometry'])
            profile.update(label=str(saved['label']), local=True,
                           instprm_filename=str(LOCAL_INSTRUMENT_DIR / f'{key}.instprm'),
                           notes='Saved on this computer.',
                           calibration_range=saved.get('calibration_range'),
                           provenance=saved.get('provenance', 'uploaded'))
            path = Path(profile['instprm_filename'])
            values = parse_instprm(path.read_bytes())
            profile.update(polariz=values.get('Polariz.', 0.5),
                           sh_l=values['SH/L'], zero_seed=values['Zero'],
                           wavelength=values.get('Lam', values.get('Lam1')),
                           spectrum='cu_doublet' if 'Lam2' in values else 'single')
            profiles[key] = profile
        except (OSError, ValueError, KeyError, TypeError) as exc:
            warnings.warn(f'Local instrument {metadata.name} could not be loaded: {exc}')
    return profiles


def validate_profile_range(values, lower, upper):
    """Reject a Gaussian variance that becomes negative inside the fit range."""
    if not (0 < lower < upper < 180):
        raise ValueError('Choose a fitting range between 0 and 180° 2θ.')
    low, high = math.tan(math.radians(lower / 2)), math.tan(math.radians(upper / 2))
    u, v, w = (values[k] for k in ('U', 'V', 'W'))
    points = [low, high]
    if u > 0 and low < -v / (2 * u) < high:
        points.append(-v / (2 * u))
    if min(u * t * t + v * t + w for t in points) < 0:
        raise ValueError('This .instprm predicts negative Gaussian variance in the fitting range. Use its calibrated range or recalibrate.')


def get_instrument_profile(key=None):
    key = key or DEFAULT_INSTRUMENT
    if key in ('generic', 'default'):
        key = DEFAULT_INSTRUMENT
    profiles = get_instrument_profiles()
    if key not in profiles:
        raise ValueError('Instrument profile is unavailable. Select a generic geometry or upload its .instprm again.')
    return profiles[key]


def instrument_file(profile):
    filename = profile.get('instprm_filename')
    return str(TOOLKIT_ROOT / filename) if filename else None


def configure_histogram_geometry(histogram, geometry):
    geometry_profile(geometry)
    sample = histogram.data['Sample Parameters']
    sample['Type'] = 'Bragg-Brentano' if geometry == 'bragg_brentano' else 'Debye-Scherrer'
    for key in ('Shift', 'DisplaceX', 'DisplaceY'):
        sample[key] = [0.0, False]
    # Single-wavelength imports default to transmission and can omit these
    # required Bragg-Brentano fields entirely.
    for key, value in (('Transparency', 0.0), ('SurfRoughA', 0.0), ('SurfRoughB', 0.0), ('Absorption', 0.0)):
        sample.setdefault(key, [value, False])


def save_local_instrument(label, geometry, content, calibration_range=None,
                          provenance='uploaded'):
    label = str(label or '').strip()
    if not label or len(label) > 100 or any(ord(c) < 32 for c in label):
        raise ValueError('Enter an instrument name of 1–100 characters.')
    geometry_profile(geometry)
    parse_instprm(content)
    if isinstance(content, bytes):
        content = content.decode('utf-8-sig')
    key = 'local_' + uuid.uuid4().hex
    LOCAL_INSTRUMENT_DIR.mkdir(parents=True, exist_ok=True)
    profile_path = LOCAL_INSTRUMENT_DIR / f'{key}.instprm'
    metadata_path = LOCAL_INSTRUMENT_DIR / f'{key}.json'
    temporary = metadata_path.with_suffix('.json.tmp')
    try:
        profile_path.write_text(content, encoding='utf-8')
        temporary.write_text(json.dumps({
            'label': label, 'geometry': geometry,
            'calibration_range': calibration_range, 'provenance': provenance,
        }, indent=2), encoding='utf-8')
        os.replace(temporary, metadata_path)
    except OSError:
        for path in (profile_path, temporary):
            path.unlink(missing_ok=True)
        raise
    return key


def silicon_640g_phase():
    """Certified Si cell and the GSAS-II Fd-3m origin setting."""
    a = 5.431109
    cif = f"""data_Si640g
_chemical_formula_sum 'Si'
_cell_formula_units_Z 8
_cell_length_a {a}
_cell_length_b {a}
_cell_length_c {a}
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_IT_number 227
_space_group_name_H-M_alt 'F d -3 m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1 Si 0.125 0.125 0.125 1
"""
    return dict(name='NIST SRM 640g Si', formula='Si', system='cubic',
                spacegroup_number=227, spacegroup='F d -3 m', a=a, b=a, c=a,
                alpha=90, beta=90, gamma=90, Z=8, cif_text=cif, source='generated')
