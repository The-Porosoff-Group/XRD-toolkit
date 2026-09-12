import copy
import tempfile
import unittest

from openpyxl import load_workbook

from modules.xrd import _write_summary_xlsx
from modules.xrd.size_reporting import apply_size_reporting, size_reporting_settings
from test_xrd_publication_exports import publication_result


def fitted_result():
    result = publication_result()
    result['phase_results'] = [result['phase_results'][0]]
    result['phase_patterns'] = [result['phase_patterns'][0]]
    result['phase_results'][0].update(
        crystallite_size_nm=8.2, crystallite_size_source='gsas_hap_size',
        gsas_hap_size_nm=8.2049, gsas_size_model='isotropic', gsas_size_lg_mix=1.0)
    return result


class SizeReportingTests(unittest.TestCase):
    def test_conversion_uses_unrounded_size_and_preserves_fit(self):
        result = fitted_result()
        original = copy.deepcopy(result)
        apply_size_reporting(result)
        phase = result['phase_results'][0]
        self.assertEqual(phase['scherrer_equivalent_size_nm'], 7.38)
        self.assertEqual(phase['gsas_hap_size_nm'], 8.2049)
        self.assertEqual(phase['crystallite_size_nm'], 8.2)
        for key in ('tt', 'y_obs', 'y_calc', 'residuals', 'statistics'):
            self.assertEqual(result[key], original[key])
        for mode, k, expected in [('hap', 1, 8.2), ('scherrer', 0.94, 7.71),
                                  ('both', 0.9, 7.38)]:
            apply_size_reporting(result, mode, k)
            self.assertEqual(phase['scherrer_equivalent_size_nm'], expected)
            self.assertEqual(phase['gsas_hap_size_nm'], 8.2049)
            self.assertEqual(len(phase['size_reporting_lines']), 2 if mode == 'both' else 1)

    def test_unverified_or_other_size_estimates_are_not_converted(self):
        for overrides in (
            {'crystallite_size_source': 'phase_pattern_scherrer'},
            {'crystallite_size_source': 'profile_Y'},
            {'crystallite_size_source': None, 'crystallite_size_nm': None},
            {'gsas_size_lg_mix': 0.5},
            {'gsas_size_lg_mix': None},
            {'gsas_size_model': 'uniaxial'},
        ):
            with self.subTest(overrides=overrides):
                result = fitted_result()
                phase = result['phase_results'][0]
                phase.update(overrides)
                apply_size_reporting(result)
                self.assertIsNone(phase['scherrer_equivalent_size_nm'])
                if phase['crystallite_size_source'] != 'gsas_hap_size':
                    self.assertIsNone(phase['gsas_hap_size_nm'])

    def test_invalid_settings_are_rejected(self):
        for k in (0, -0.9, float('nan'), float('inf'), '', 'abc', None):
            with self.subTest(k=k), self.assertRaises(ValueError):
                size_reporting_settings('both', k)
        with self.assertRaises(ValueError):
            size_reporting_settings('unknown', 0.9)

    def test_workbook_obeys_selected_convention_and_records_method(self):
        with tempfile.TemporaryDirectory() as directory:
            for mode in ('hap', 'scherrer', 'both'):
                with self.subTest(mode=mode):
                    result = apply_size_reporting(fitted_result(), mode, 0.9)
                    path = _write_summary_xlsx(result, {'sample_id': mode},
                                               'GSAS-II', directory)
                    workbook = load_workbook(path, read_only=True, data_only=True)
                    try:
                        rows = {row[0]: row[1] for row in workbook['Summary'].values}
                        hap_key = 'GSAS HAP size (nm; K=1)'
                        scherrer_key = 'Scherrer-equivalent size (nm)'
                        self.assertEqual(hap_key in rows, mode in ('hap', 'both'))
                        self.assertEqual(scherrer_key in rows, mode in ('scherrer', 'both'))
                        if hap_key in rows:
                            self.assertAlmostEqual(rows[hap_key], 8.2049)
                        if scherrer_key in rows:
                            self.assertEqual(rows[scherrer_key], 7.38)
                            self.assertEqual(rows['Scherrer shape factor K'], 0.9)
                        self.assertIn('not an independent', rows['Size reporting method'])
                    finally:
                        workbook.close()


if __name__ == '__main__':
    unittest.main()
