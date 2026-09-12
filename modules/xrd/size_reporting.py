"""Report size conventions without modifying any fitted parameter or curve."""
import math


def size_reporting_settings(mode='both', scherrer_k=0.9):
    if mode not in ('hap', 'scherrer', 'both'):
        raise ValueError('Size reporting must be GSAS HAP, Scherrer-equivalent, or Both.')
    try:
        k = float(scherrer_k)
    except (ValueError, TypeError):
        raise ValueError('Scherrer K must be a finite positive number.') from None
    if not math.isfinite(k) or k <= 0:
        raise ValueError('Scherrer K must be a finite positive number.')
    return {'mode': mode, 'scherrer_k': k}


def apply_size_reporting(result, mode='both', scherrer_k=0.9):
    """Attach explicit GSAS and K-adjusted sizes; retain native size fields.

    GSAS-II's isotropic, Lorentzian HAP Size uses Scherrer K=1. Rescaling
    that fitted size by K is not an independent single-reflection fit.
    Other profile types and fallback estimates must not be relabeled HAP.
    """
    settings = size_reporting_settings(mode, scherrer_k)
    k = settings['scherrer_k']
    result['size_reporting'] = settings
    for phase in result.get('phase_results', []):
        is_hap = phase.get('crystallite_size_source') == 'gsas_hap_size'
        value = phase.get('gsas_hap_size_nm')
        if value is None:
            value = phase.get('crystallite_size_nm')
        try:
            hap = float(value) if is_hap else None
            if hap is not None and (not math.isfinite(hap) or hap <= 0):
                hap = None
        except (TypeError, ValueError):
            hap = None
        phase['gsas_hap_size_nm'] = hap
        phase['scherrer_k'] = k
        phase['size_reporting_mode'] = mode
        equivalent = None
        if hap is None:
            note = ('No refined GSAS HAP size. Enable Size in the phase card '
                    'and rerun GSAS-II to obtain a Scherrer-equivalent value.')
        elif (phase.get('gsas_size_model') != 'isotropic'
              or phase.get('gsas_size_lg_mix') != 1.0):
            note = ('Scherrer-equivalent conversion requires a verified isotropic '
                    'Lorentzian size model. Rerun the fit if model metadata is missing.')
        else:
            if not math.isfinite(k * hap):
                raise ValueError('Scherrer K is too large for this fitted size.')
            equivalent = round(k * hap, 2)
            note = (f'Scherrer-equivalent = {k:g} × GSAS HAP size (K=1), '
                    'using the fitted size broadening. This is not an independent '
                    'single-peak Scherrer measurement. Instrument calibration and '
                    'size/strain model assumptions still apply.')
        phase['scherrer_equivalent_size_nm'] = equivalent
        phase['size_reporting_note'] = note
        lines = []
        if mode in ('hap', 'both'):
            lines.append({'label': 'GSAS HAP size (K = 1)',
                          'value_nm': round(hap, 2) if hap is not None else None})
        if mode in ('scherrer', 'both'):
            lines.append({'label': f'Scherrer-equivalent size (K = {k:g})',
                          'value_nm': equivalent})
        if hap is None and phase.get('crystallite_size_nm') is not None:
            lines.append({'label': 'Other size estimate (unconverted)',
                          'value_nm': phase['crystallite_size_nm']})
        phase['size_reporting_lines'] = lines
    return result
