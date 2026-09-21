"""Numerical checks for mass quantification independent of display intensity."""

import copy
import math
import unittest

from modules.xrd.quantification import (
    attach_phase_area_diagnostics,
    official_mass_fraction_uncertainties,
    phase_area_diagnostics,
    resolve_mass_fractions,
)


class XrdMassFractionTests(unittest.TestCase):
    names = ['light_cell', 'heavy_cell']

    def resolve(self, scales=None, masses=None, official=None):
        return resolve_mass_fractions(
            self.names,
            scales if scales is not None else dict.fromkeys(self.names, 1.0),
            masses if masses is not None else dict(zip(self.names, [100., 300.])),
            official_mass_fractions=official,
        )

    def test_equal_cell_counts_are_not_equal_masses(self):
        result = self.resolve()
        self.assertAlmostEqual(result['fractions_pct']['light_cell'], 25.)
        self.assertAlmostEqual(result['fractions_pct']['heavy_cell'], 75.)

    def test_scale_and_mass_compensate_for_equal_sample_masses(self):
        result = self.resolve(scales={'light_cell': 3., 'heavy_cell': 1.})
        self.assertAlmostEqual(result['fractions_pct']['light_cell'], 50.)
        self.assertAlmostEqual(result['fractions_pct']['heavy_cell'], 50.)

    def test_overall_intensity_scale_does_not_change_composition(self):
        original = self.resolve(scales={'light_cell': 3., 'heavy_cell': 2.})
        rescaled = self.resolve(scales={'light_cell': 30., 'heavy_cell': 20.})
        self.assertEqual(original['fractions_pct'], rescaled['fractions_pct'])

    def test_missing_or_invalid_mass_never_becomes_raw_scale_percentage(self):
        for mass in (None, 0., -1., math.nan, math.inf):
            with self.subTest(mass=mass):
                result = self.resolve(masses={'light_cell': mass, 'heavy_cell': 300.})
                self.assertTrue(all(v is None for v in result['fractions_pct'].values()))
                self.assertTrue(result['note'])

    def test_invalid_scales_or_zero_total_do_not_create_composition(self):
        for scales in ({'light_cell': -1., 'heavy_cell': 2.},
                       {'light_cell': math.nan, 'heavy_cell': 1.},
                       {'light_cell': math.inf, 'heavy_cell': 1.},
                       {'light_cell': 0., 'heavy_cell': 0.}):
            with self.subTest(scales=scales):
                result = self.resolve(scales=scales)
                self.assertTrue(all(v is None for v in result['fractions_pct'].values()))

    def test_zero_abundance_of_one_phase_is_valid(self):
        result = self.resolve(scales={'light_cell': 0., 'heavy_cell': 1.})
        self.assertEqual(result['fractions_pct'], {'light_cell': 0., 'heavy_cell': 100.})

    def test_complete_gsas_result_is_authoritative_and_converts_fraction_units(self):
        result = self.resolve(official={'light_cell': (.25, .03), 'heavy_cell': (.75, .03)})
        self.assertEqual(result['fractions_pct'], {'light_cell': 25., 'heavy_cell': 75.})
        self.assertEqual(result['sigmas_pct'], {'light_cell': 3., 'heavy_cell': 3.})
        self.assertEqual(result['method'], 'gsasii_compute_mass_fracs')

    def test_complete_official_result_can_supply_unreadable_mass_metadata(self):
        result = self.resolve(masses={}, official={'light_cell': (.4, .03), 'heavy_cell': (.6, .03)})
        self.assertEqual(result['fractions_pct'], {'light_cell': 40., 'heavy_cell': 60.})

    def test_incomplete_or_invalid_official_result_is_not_mixed_with_fallback(self):
        for official in ('not a mass result', ['invalid'], {'light_cell': (.9, .01)},
                         {'light_cell': (.4, .01), 'heavy_cell': (.6, .01)},
                         {'light_cell': (.9, .01), 'heavy_cell': (.9, .01)},
                         {'light_cell': (math.nan, .01), 'heavy_cell': (.75, .01)},
                         {'light_cell': (-.1, .01), 'heavy_cell': (1.1, .01)}):
            with self.subTest(official=official):
                result = self.resolve(official=official)
                self.assertEqual(result['fractions_pct'], {'light_cell': 25., 'heavy_cell': 75.})
                self.assertEqual(result['sigmas_pct'], {})

    def test_area_integrates_over_angle_even_on_a_nonuniform_grid(self):
        # Areas are 4 and 1, although the unweighted intensity sums are 3 and 2.
        diagnostics = phase_area_diagnostics([0., 1., 4.], [[1., 1., 1.], [2., 0., 0.]], 'gsasii_isolation')
        self.assertAlmostEqual(diagnostics[0]['integrated_phase_fraction_%'], 80.)
        self.assertAlmostEqual(diagnostics[1]['integrated_phase_fraction_%'], 20.)
        for phase in diagnostics:
            self.assertFalse(any(key.startswith('weight_fraction') for key in phase))

    def test_area_diagnostic_does_not_mutate_input_patterns(self):
        patterns = [[1., 2., 3.], [3., 2., 1.]]
        original = copy.deepcopy(patterns)
        phase_area_diagnostics([0., 1., 2.], patterns, 'profile_reconstruction')
        self.assertEqual(patterns, original)

    def test_drawing_methods_cannot_overwrite_mass_or_its_uncertainty(self):
        for method, patterns in (
                ('gsasii_isolation', [[4., 4., 4.], [1., 1., 1.]]),
                ('profile_reconstruction', [[1., 1., 1.], [4., 4., 4.]]),
                ('equal_split', [[1., 1., 1.], [1., 1., 1.]])):
            with self.subTest(method=method):
                phases = [
                    {'weight_fraction_%': 25., 'weight_fraction_err_%': 3.,
                     'weight_fraction_sigma_propagated_%': 3.,
                     'weight_fraction_method': 'gsasii_compute_mass_fracs',
                     'weight_fraction_err_source': 'gsasii_calcMassFracs'},
                    {'weight_fraction_%': 75., 'weight_fraction_err_%': 3.,
                     'weight_fraction_sigma_propagated_%': 3.,
                     'weight_fraction_method': 'gsasii_compute_mass_fracs',
                     'weight_fraction_err_source': 'gsasii_calcMassFracs'},
                ]
                mass_results = copy.deepcopy(phases)
                attach_phase_area_diagnostics(phases, [0., 1., 2.], patterns, method)
                for phase, mass_result in zip(phases, mass_results):
                    for key, value in mass_result.items():
                        self.assertEqual(phase[key], value)
                    self.assertEqual(phase['integrated_phase_fraction_method'], method)
                    self.assertNotEqual(phase['integrated_phase_fraction_%'], phase['weight_fraction_%'])

    def test_zero_signal_and_invalid_grids_have_no_area_percentage(self):
        for angles, patterns in (([0., 1.], [[0., 0.], [0., 0.]]),
                                 ([0., 0.], [[1., 1.], [1., 1.]]),
                                 ([1., 0.], [[1., 1.], [1., 1.]]),
                                 ([0., 1.], [[math.nan, 1.], [1., 1.]])):
            with self.subTest(angles=angles, patterns=patterns):
                diagnostics = phase_area_diagnostics(angles, patterns, 'test')
                self.assertTrue(all(d['integrated_phase_fraction_%'] is None for d in diagnostics))

    def test_official_uncertainties_work_without_individual_scale_esds(self):
        uncertainty = official_mass_fraction_uncertainties(
            {'light_cell': 25., 'heavy_cell': 75.},
            {'light_cell': 3., 'heavy_cell': 3.})
        for phase in uncertainty.values():
            self.assertEqual(phase['weight_fraction_err_%'], 3.)
            self.assertEqual(phase['weight_fraction_sigma_propagated_%'], 3.)
            self.assertIn('gsasii_calcMassFracs', phase['weight_fraction_err_source'])

    def test_reporting_floor_preserves_formal_sigma_and_binary_closure(self):
        uncertainty = official_mass_fraction_uncertainties(
            {'light_cell': 25., 'heavy_cell': 75.},
            {'light_cell': .05, 'heavy_cell': .05})
        for phase in uncertainty.values():
            self.assertEqual(phase['weight_fraction_err_%'], 1.5)
            self.assertEqual(phase['weight_fraction_sigma_propagated_%'], .05)
            self.assertIn('systematic_floor', phase['weight_fraction_err_source'])


if __name__ == '__main__':
    unittest.main()
