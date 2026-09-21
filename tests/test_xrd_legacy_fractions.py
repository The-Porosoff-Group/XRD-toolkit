"""Legacy fitted intensity areas must never be reported as weight percent."""

import copy
import unittest

import numpy as np

from modules.xrd.lebail import (
    _integrated_phase_fractions,
    run_lebail,
    run_rietveld,
)


class LegacyFractionTests(unittest.TestCase):
    @staticmethod
    def _phases():
        return [
            {'name': 'Ni', 'formula': 'Ni', 'Z': 1, 'a': 3.5,
             'system': 'cubic', 'spacegroup_number': 1,
             'sites': [('Ni', 0.0, 0.0, 0.0, 1.0)]},
            {'name': 'Cu', 'formula': 'Cu', 'Z': 1, 'a': 4.1,
             'system': 'cubic', 'spacegroup_number': 1,
             'sites': [('Cu', 0.0, 0.0, 0.0, 1.0)]},
        ]

    def test_both_legacy_engines_withhold_mass_for_complete_and_missing_metadata(self):
        tt = np.linspace(20.0, 80.0, 601)
        observed = 50.0 + 1000.0 * np.exp(-((tt - 44.0) / 0.4) ** 2)
        for engine, iteration_option, reason in (
            (run_lebail, 'max_outer', 'unavailable_lebail_free_intensities'),
            (run_rietveld, 'max_iter', 'unavailable_legacy_profile_normalization'),
        ):
            baseline_areas = None
            for missing in (None, 'Z', 'formula'):
                with self.subTest(engine=engine.__name__, missing=missing):
                    phases = copy.deepcopy(self._phases())
                    if missing:
                        phases[0].pop(missing)
                    # No optimizer is needed to exercise the final reporting
                    # path; keep these tests deterministic and inexpensive.
                    result = engine(
                        tt, observed, None, phases, 1.54056,
                        **{iteration_option: 0})
                    fitted_phases = result['phase_results']
                    areas = [phase['integrated_phase_fraction_%']
                             for phase in fitted_phases]
                    self.assertAlmostEqual(sum(areas), 100.0, places=1)
                    self.assertTrue(all(area > 0 for area in areas))
                    if baseline_areas is None:
                        baseline_areas = areas
                    else:
                        self.assertEqual(areas, baseline_areas)
                    for phase in fitted_phases:
                        self.assertIsNone(phase['weight_fraction_%'])
                        self.assertEqual(phase['weight_fraction_method'], reason)
                        self.assertIn('GSAS-II', phase['weight_fraction_note'])
                        self.assertEqual(
                            phase['integrated_phase_fraction_method'],
                            'integrated_fitted_phase_intensity_area')

    def test_area_diagnostic_accounts_for_nonuniform_angle_spacing(self):
        # The first triangle has area 1.5; the second has area 0.5.
        # Merely summing the sampled intensities would incorrectly give 50/50.
        fractions = _integrated_phase_fractions(
            [0.0, 1.0, 3.0], [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
        self.assertEqual(fractions, [75.0, 25.0])

    def test_zero_total_intensity_has_no_fraction(self):
        self.assertEqual(
            _integrated_phase_fractions([0.0, 1.0], [[0.0, 0.0], [0.0, 0.0]]),
            [None, None])


if __name__ == '__main__':
    unittest.main()
