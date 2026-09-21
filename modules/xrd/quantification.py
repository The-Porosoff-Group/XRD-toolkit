"""Keep crystalline mass fractions separate from diffraction-area diagnostics."""

import math
from collections.abc import Mapping

import numpy as np


def _finite_number(value, minimum=0.0, strictly_positive=False):
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    if not math.isfinite(value) or value < minimum:
        return None
    if strictly_positive and value <= 0.0:
        return None
    return value


def resolve_mass_fractions(phase_names, scales, unit_cell_masses,
                           official_mass_fractions=None):
    """Return one coherent all-phase mass result, never normalized raw scales.

    GSAS-II HAP scales multiply unit-cell mass (no extra cell-volume factor).
    ``official_mass_fractions`` is ComputeMassFracs' mapping of phase names to
    (mass fraction, sigma), both in fractional units. A complete normalized
    official result may be used when independently readable masses are absent.
    Invalid/missing scales always make quantification unavailable.
    """
    names = list(phase_names)
    unavailable = {
        'fractions_pct': dict.fromkeys(names), 'sigmas_pct': {},
        'mass_weighted_values': {}, 'method': 'unavailable', 'note': None,
    }
    clean_scales = {name: _finite_number(scales.get(name)) for name in names}
    if (not names or len(set(names)) != len(names)
            or any(value is None for value in clean_scales.values())
            or not any(value > 0.0 for value in clean_scales.values())):
        unavailable['note'] = (
            'Weight fractions unavailable: phase scales must be finite, '
            'nonnegative, and have a positive total.')
        return unavailable

    clean_masses = {
        name: _finite_number(unit_cell_masses.get(name), strictly_positive=True)
        for name in names
    }
    weighted = {}
    calculated = {}
    if all(value is not None for value in clean_masses.values()):
        weighted = {name: clean_scales[name] * clean_masses[name]
                    for name in names}
        total = sum(weighted.values())
        if math.isfinite(total) and total > 0.0:
            calculated = {name: weighted[name] / total * 100.0 for name in names}
        else:
            weighted = {}

    official = {}
    sigmas = {}
    supplied = (official_mass_fractions
                if isinstance(official_mass_fractions, Mapping) else {})
    for name in names:
        pair = supplied.get(name)
        if not isinstance(pair, (list, tuple)) or not pair:
            break
        fraction = _finite_number(pair[0])
        if fraction is None or fraction > 1.0:
            break
        official[name] = fraction * 100.0
        sigma = _finite_number(pair[1]) if len(pair) > 1 else None
        if sigma is not None:
            sigmas[name] = sigma * 100.0
    coherent_official = (
        len(official) == len(names)
        and math.isclose(sum(official.values()), 100.0, rel_tol=0.0, abs_tol=1e-5)
        and all(clean_scales[name] != 0.0 or official[name] == 0.0
                for name in names)
        and (not calculated or all(math.isclose(
            official[name], calculated[name], rel_tol=1e-5, abs_tol=1e-5)
            for name in names))
    )
    if coherent_official:
        return {
            'fractions_pct': official, 'sigmas_pct': sigmas,
            'mass_weighted_values': weighted,
            'method': 'gsasii_compute_mass_fracs',
            'note': 'Mass percentages normalized over the modeled crystalline phases.',
        }
    if calculated:
        return {
            'fractions_pct': calculated, 'sigmas_pct': {},
            'mass_weighted_values': weighted, 'method': 'gsasii_mass_fraction',
            'note': 'Mass percentages from GSAS-II Scale × unit-cell mass, '
                    'normalized over the modeled crystalline phases.',
        }
    unavailable['note'] = (
        'Weight fractions unavailable: valid unit-cell masses are required for '
        'every phase, or GSAS-II must provide a complete, normalized mass result.')
    return unavailable


def official_mass_fraction_uncertainties(fractions_pct, sigmas_pct):
    """Apply existing reporting floors to official sigmas without scale ESDs.

    All values are percentage points; the formal propagated sigma remains
    separately available from the heuristic systematic reporting floor.
    """
    result = {}
    for name, fraction in fractions_pct.items():
        fraction = _finite_number(fraction)
        sigma = _finite_number(sigmas_pct.get(name))
        if fraction is None or sigma is None:
            continue
        floor = max(1.0, 0.02 * fraction)
        reported = max(sigma, floor)
        source = 'gsasii_calcMassFracs'
        if reported > sigma * 1.01:
            source += '+systematic_floor'
        result[name] = {
            'weight_fraction_err_%': reported,
            'weight_fraction_sigma_propagated_%': sigma,
            'weight_fraction_systematic_floor_%': floor,
            'weight_fraction_err_source': source,
        }
    if len(fractions_pct) == 2 and len(result) == 2:
        binary_sigma = max(item['weight_fraction_err_%'] for item in result.values())
        for item in result.values():
            item['weight_fraction_err_%'] = binary_sigma
            item['weight_fraction_err_source'] += '+binary_closure'
    return result


def phase_area_diagnostics(two_theta, phase_patterns, method):
    """Integrate displayed phase intensity over angle; do not return mass keys.

    Trapezoidal integration handles uneven angular spacing. Reconstructed or
    equal-split patterns remain display diagnostics, regardless of their area.
    """
    result = [{
        'integrated_phase_fraction_%': None,
        'integrated_phase_fraction_method': method,
        'integrated_phase_fraction_note': (
            'Diffraction-area share over the fitted angular range; not a mass percentage.'),
    } for _ in phase_patterns]
    try:
        angles = np.asarray(two_theta, dtype=float)
        patterns = [np.asarray(pattern, dtype=float) for pattern in phase_patterns]
        if (angles.ndim != 1 or len(angles) < 2
                or not np.all(np.isfinite(angles))
                or not np.all(np.diff(angles) > 0.0)
                or any(pattern.shape != angles.shape
                       or not np.all(np.isfinite(pattern))
                       or np.any(pattern < 0.0) for pattern in patterns)):
            raise ValueError('invalid angular grid or phase intensity')
        integrate = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
        areas = [float(integrate(pattern, angles)) for pattern in patterns]
        total = sum(areas)
        if not math.isfinite(total) or total <= 0.0:
            raise ValueError('no positive finite integrated phase intensity')
        for item, area in zip(result, areas):
            item['integrated_phase_fraction_%'] = round(area / total * 100.0, 1)
    except (TypeError, ValueError, OverflowError):
        for item in result:
            item['integrated_phase_fraction_note'] = (
                'Diffraction-area share unavailable: invalid or zero integrated intensity.')
    return result


def attach_phase_area_diagnostics(phase_results, two_theta, phase_patterns, method):
    """Attach display-area diagnostics without replacing any mass result or ESD."""
    for phase, diagnostic in zip(
            phase_results, phase_area_diagnostics(two_theta, phase_patterns, method)):
        phase.update(diagnostic)
        area_pct = diagnostic['integrated_phase_fraction_%']
        mass_pct = _finite_number(phase.get('weight_fraction_%'))
        difference = (round(area_pct - mass_pct, 2)
                      if area_pct is not None and mass_pct is not None else None)
        phase['integrated_minus_weight_fraction_pp'] = difference
        phase['integrated_minus_hh_fraction_pp'] = difference
    return phase_results
