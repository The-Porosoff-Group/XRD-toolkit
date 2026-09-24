"""Readable interface settings for repeating an XRD fit."""

import os

import numpy as np


def json_value(value):
    """Copy fit inputs into plain JSON values without rounding numbers."""
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    return value


_MISSING = 'Not recorded'
_GLOBAL_CONTROLS = (
    ('xrd-verification-mode', 'verification_mode', 'Quick constrained fit'),
    ('xrd-verify-cell', 'verify_refine_cell', 'Refine cell too'),
    ('xrd-phase-isolation', 'phase_isolation', 'Phase isolation'),
    ('xrd-zero-not-disp', 'verify_use_zero_not_displace', 'Refine Zero (fix Disp)'),
    ('xrd-refine-x', 'verify_refine_x', 'Free X'),
    ('xrd-fix-y', 'verify_fix_y', 'Fix Y'),
    ('xrd-y-nonneg', 'verify_y_nonnegative', 'Y ≥ 0'),
    ('xrd-refine-uiso', 'verify_refine_uiso', 'Refine Uiso'),
    ('xrd-refine-xyz', 'refine_xyz', 'Refine XYZ'),
)


def _number(value):
    if value is None or value == '':
        return _MISSING
    try:
        number = float(value)
        return number if np.isfinite(number) else _MISSING
    except (TypeError, ValueError):
        return _MISSING


def _checked(value):
    if isinstance(value, (bool, np.bool_)):
        return 'Checked' if value else 'Unchecked'
    return _MISSING


def _axis(value):
    if isinstance(value, (list, tuple)):
        return ' '.join(str(item) for item in value)
    return str(value) if value not in (None, '') else 'Blank (automatic)'


def _options(value, index):
    if isinstance(value, list):
        return value[index] if index < len(value) and isinstance(value[index], dict) else {}
    if isinstance(value, dict):
        item = value.get(str(index), value.get(index, {}))
        return item if isinstance(item, dict) else {}
    return {}


def fit_parameter_rows(result, metadata, method_label):
    """Export a small, explicit list of initial inputs and visible UI settings.

    Never flatten the native project, arrays, or CIF text into the workbook.
    Missing older settings remain unknown instead of borrowing current defaults
    or presenting refined results as the original phase-card inputs.
    """
    settings = result.get('fit_settings') or {}
    submitted = settings.get('submitted_parameters') or {}
    effective = settings.get('effective_parameters') or {}
    interface = settings.get('interface_settings') or {}
    controls = interface.get('controls') or {}
    checks = controls.get('checkboxes') or {}
    native = result.get('gsas_native_parameters') or {}
    resolved = native.get('effective_settings') or {}
    rows = []

    def add(section, label, value):
        if value is None or value == '':
            value = _MISSING
        rows.append({'Section': section, 'Setting': label, 'Value': value})

    def entered(key, submitted_key=None):
        return controls.get(key, submitted.get(submitted_key or key))

    add('Repeat this fit', 'Instructions',
        'Use the same scan, phase cards/CIF files, instrument file, and settings below. '
        'Use matching toolkit and GSAS-II versions.')
    add('Repeat this fit', 'GSAS-II project', settings.get('project_file') or
        ('See the companion GPX for the fitted GSAS-II model.' if native.get('native_project')
         else 'No companion GPX recorded for this result.'))
    add('Run', 'Toolkit version', (settings.get('software') or {}).get('toolkit_version'))
    version = (native.get('software') or {}).get('gsas_version')
    if isinstance(version, dict):
        version = version.get('git_versiontag') or version.get('git_version')
    if version:
        add('Run', 'GSAS-II version', str(version))
    add('Run', 'Method', method_label)
    add('Run', 'Data file', (settings.get('source') or {}).get('filename') or metadata.get('source_file'))
    add('Run', 'Sample ID', metadata.get('sample_id'))
    if not submitted and not controls:
        add('Repeat this fit', 'Settings unavailable',
            'These interface settings were not captured in this older result. Run the fit again to record them.')
        return rows

    instrument = entered('instrument') or effective.get('instrument')
    labels = {'smartlab':'Rigaku SmartLab (BB)', 'synergy_s':'Synergy-S (capillary)',
              'benchtop_cu':'Benchtop Cu — flat plate (Si 640g)', 'none':'None / calibration',
              'generic_flat_plate':'Generic flat plate', 'generic_capillary':'Generic capillary',
              'auto':'Auto (legacy selection)'}
    add('Scan settings', 'Instrument', labels.get(instrument, instrument))
    if instrument == 'auto':
        actual = effective.get('instrument')
        add('Scan settings', 'Instrument used', labels.get(actual, actual))
    profile = settings.get('instrument_file') or {}
    add('Scan settings', 'Instrument file', profile.get('filename') or
        (os.path.basename(submitted['instprm_file']) if submitted.get('instprm_file')
         else (resolved.get('instrument_profile') or {}).get('instprm_filename')))
    if profile.get('sha256'):
        add('Scan settings', 'Instrument file SHA-256', profile['sha256'])
    add('Scan settings', 'X-ray Source', interface.get('wavelength_label') or submitted.get('wavelength_label'))
    add('Scan settings', 'Wavelength entered (Å)', _number(entered('wavelength')))
    actual_wavelength = effective.get('wavelength')
    if actual_wavelength is not None:
        add('Scan settings', 'Wavelength used (Å)', _number(actual_wavelength))
    add('Scan settings', '2θ Min (°)', _number(entered('tt_min')))
    add('Scan settings', '2θ Max (°)', _number(entered('tt_max')))
    background = entered('n_bg_coeffs')
    add('Scan settings', 'Background terms', 'Auto' if str(background).lower() == 'auto' else _number(background))
    if resolved.get('n_bg_coeffs') is not None:
        add('Scan settings', 'Background terms used', _number(resolved['n_bg_coeffs']))

    if 'gsas' in method_label.lower():
        for dom_id, key, label in _GLOBAL_CONTROLS:
            add('GSAS controls', label, _checked(checks.get(dom_id, submitted.get(key))))
        fixed_y = entered('fix_y_value', 'verify_y_fixed_value')
        blank_y_recorded = ('fix_y_value' in controls or
                            submitted.get('verify_fix_y') is True)
        add('GSAS controls', 'Fix Y value (centideg)',
            ('Blank (use instrument file)' if blank_y_recorded else _MISSING)
            if fixed_y in (None, '') else _number(fixed_y))

    phases = settings.get('input_phases') or submitted.get('phases') or []
    for index, phase in enumerate(phases):
        section = f"Phase {index + 1}: {phase.get('name') or phase.get('formula') or 'unnamed'}"
        add(section, 'Formula', phase.get('formula'))
        source = phase.get('source')
        add(section, 'Phase source', {'mp':'Materials Project', 'cod':'Crystallography Open Database',
                                    'manual':'Uploaded CIF', 'generated':'Generated carbide CIF'}.get(source, source))
        generated = source == 'generated' or bool(phase.get('generated_cif_model'))
        card_id = (phase.get('cod_id') if generated else
                   phase.get('mp_id') or phase.get('cod_id'))
        add(section, 'Phase card ID', 'Manual CIF' if card_id == 'manual' else card_id)
        if generated:
            add(section, 'Reference Materials Project card', phase.get('mp_id'))
            add(section, 'Carbon vacancy fraction x', _number(phase.get('gamma_vacancy_x')))
            add(section, 'Fixed carbon occupancy', _number(phase.get('gamma_c_occupancy')))
        if phase.get('cif_filename'):
            add(section, 'CIF file', phase['cif_filename'])
        add(section, 'Space group', phase.get('spacegroup') or phase.get('spacegroup_number'))
        add(section, 'Crystal system', phase.get('system'))
        add(section, 'Z (formula units per cell)', _number(phase.get('Z')))
        for key in ('a', 'b', 'c'):
            add(section, f'Starting {key} (Å)', _number(phase.get(key)))
        angles = [phase.get(key) for key in ('alpha', 'beta', 'gamma')]
        if all(value is not None for value in angles):
            add(section, 'Starting α, β, γ (°)', ', '.join(str(value) for value in angles))
        options = dict(_options(submitted.get('phase_options'), index))
        ui_options = _options(interface.get('phase_options'), index)
        options.update(ui_options)
        for key, label in [('refine_cell','Cell'), ('refine_size','Size'), ('refine_mustrain','Mustrain')]:
            add(section, label, _checked(options.get(key)))
        add(section, 'PO', options.get('po_mode'))
        add(section, 'PO axis', _axis(ui_options.get('po_axis', options.get('po_axis_input', options.get('po_axis'))))
            if 'po_axis' in options or 'po_axis_input' in options else _MISSING)
        add(section, 'PO value', _number(options.get('po_value')))
        formula = ''.join(str(phase.get('formula') or phase.get('name') or '').lower().split())
        if 'w2c' in formula or formula == 'cw2' or 'c1w2' in formula:
            add(section, 'Uniform cell', _checked(options.get('uniform_cell')))

    add('Figure and reporting', 'Size reporting', {'both':'Both', 'hap':'GSAS HAP size (K = 1)',
        'scherrer':'Scherrer-equivalent size'}.get(entered('size_reporting_mode'), entered('size_reporting_mode')))
    add('Figure and reporting', 'Scherrer shape factor K', _number(entered('scherrer_k')))
    add('Figure and reporting', 'Figure title', metadata.get('figure_title') or 'Blank (use Sample ID)')
    add('Figure and reporting', 'Show figure title', _checked(metadata.get('show_figure_title')))
    add('Figure and reporting', 'Legend location at fit completion', entered('legend_location') or metadata.get('legend_location'))
    return rows
