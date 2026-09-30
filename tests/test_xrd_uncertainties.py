"""Regression checks for measured uncertainties and native GSAS statistics."""
import math
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np

from modules.xrd import parse_xrd_file
from modules.xrd.gsasii_backend import (
    _gsas_fit_statistics, _uncertainty_warning, _write_xye,
)
from modules.xrd.instrument_profiles import get_instrument_profiles


class UncertaintyTests(unittest.TestCase):
    def parse(self, text, suffix='.csv'):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / ('scan' + suffix)
            path.write_text(text, encoding='utf-8')
            return parse_xrd_file(path)

    def test_synergy_csv_uses_sigx_not_d_spacing(self):
        parsed = self.parse('2thetadeg,intx,d-value,sigx,count\n'
                            '0.001,0,88267.7,0,0\n'
                            '20,100,4.4,0.2,123\n30,144,2.9,0.3,456\n')
        np.testing.assert_allclose(parsed['sigma'], [.2, .3])
        np.testing.assert_allclose(parsed['intensity'], [100, 144])
        self.assertEqual(parsed['metadata']['columns']['sigma'], 'sigx')
        self.assertEqual(parsed['metadata']['excluded_unmeasured_rows'], 1)

    def test_named_columns_can_be_reordered_and_quoted(self):
        parsed = self.parse('count,SIGX,"2theta (deg)",d-value,INTX\n'
                            '99,0.7,20,4.4,50\n100,0.8,21,4.2,60\n')
        np.testing.assert_allclose(parsed['tt'], [20, 21])
        np.testing.assert_allclose(parsed['intensity'], [50, 60])
        np.testing.assert_allclose(parsed['sigma'], [.7, .8])

    def test_d_spacing_is_not_sigma_when_error_column_is_absent(self):
        parsed = self.parse('2theta,intensity,d-value\n20,100,4.4\n30,144,2.9\n')
        np.testing.assert_allclose(parsed['sigma'], [10, 12])
        self.assertEqual(parsed['metadata']['sigma_source'], 'sqrt_intensity_estimate')

    def test_zero_counts_are_valid_when_count_is_the_intensity(self):
        parsed = self.parse('2theta,counts\n20,0\n30,4\n')
        np.testing.assert_allclose(parsed['intensity'], [0, 4])
        np.testing.assert_allclose(parsed['sigma'], [1, 2])

    def test_invalid_supplied_sigma_is_not_silently_replaced(self):
        for sigma in ('0', '-1', 'nan', 'inf', ''):
            with self.subTest(sigma=sigma), self.assertRaisesRegex(ValueError, 'uncertainty'):
                self.parse(f'2theta,intensity,sigx,count\n20,10,{sigma},5\n')

    def test_headerless_xy_and_xye_remain_compatible(self):
        for suffix, text, sigma in [('.xy', '20 100\n30 144\n', [10,12]),
                                    ('.xye', '20 100 2\n30 144 3\n', [2,3]),
                                    ('.csv', '20,100,2\n30,144,3\n', [2,3])]:
            with self.subTest(suffix=suffix):
                np.testing.assert_allclose(self.parse(text,suffix)['sigma'],sigma)

    def test_no_instrument_preset_rescales_uncertainties(self):
        for key, profile in get_instrument_profiles().items():
            with self.subTest(instrument=key):
                self.assertEqual(profile['sigma_inflation_K'], 1.0)

    def test_explicit_sigx_does_not_trigger_estimated_uncertainty_warning(self):
        parsed = self.parse('2thetadeg,intx,d-value,sigx,count\n'
                            '20,10.25,4.4,0.2,123\n30,14.75,2.9,0.3,456\n')
        self.assertIsNone(_uncertainty_warning(
            parsed['intensity'], parsed['sigma'], parsed['metadata']['sigma_source']))

    def test_named_csv_sqrt_estimate_warns_for_fractional_intensity(self):
        parsed = self.parse('2theta,intensity\n20,10.25\n30,14.75\n')
        warning = _uncertainty_warning(
            parsed['intensity'], parsed['sigma'], parsed['metadata']['sigma_source'])
        self.assertIn('estimated as sqrt(max(I, 1))', warning)

    def test_generic_xy_and_xye_record_uncertainty_source(self):
        for suffix, text, expected, warns in (
                ('.xy', '20 10.25\n30 14.75\n', 'sqrt_intensity_estimate', True),
                ('.xy', '20 10\n30 14\n', 'sqrt_intensity_estimate', False),
                ('.xye', '20 10.25 0.2\n30 14.75 0.3\n', 'input_column', False)):
            with self.subTest(suffix=suffix, text=text):
                parsed = self.parse(text, suffix)
                self.assertEqual(parsed['metadata']['sigma_source'], expected)
                self.assertEqual(bool(_uncertainty_warning(
                    parsed['intensity'], parsed['sigma'], expected)), warns)

    def test_direct_backend_missing_sigma_warns_but_supplied_sigma_does_not(self):
        self.assertIsNotNone(_uncertainty_warning([10.25, 14.75], None))
        self.assertIsNone(_uncertainty_warning([10.25, 14.75], [0.2, 0.3]))

    def test_gsas_xye_preserves_small_positive_uncertainties(self):
        values = np.array([[20.123456789012345, 0.0000123456789012345,
                            0.000000123456789012345],
                           [30.987654321098765, 12345.678901234567,
                            1.2345678901234567]])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'scan.xye'
            _write_xye(path, values[:, 0], values[:, 1], values[:, 2])
            written = np.loadtxt(path)
        np.testing.assert_array_equal(written, values)
        self.assertTrue(np.all(written[:, 2] > 0))
        np.testing.assert_array_equal(1 / written[:, 2] ** 2,
                                      1 / values[:, 2] ** 2)

    def histogram(self, factor=1.0):
        # Two fitted residuals: 2/2 and -3/1. A masked large residual and
        # another outside the fit range must not contribute to the sum.
        arrays = dict(x=np.array([20.,30.,40.,50.]),
                      yobs=np.ma.array([10.,20.,999.,999.], mask=[0,0,1,0]),
                      ycalc=np.array([8.,23.,0.,0.]),
                      yweight=np.array([.25,1.,1.,1.]))
        return SimpleNamespace(data={'data':[{'wtFactor':factor}],
                                     'Limits':[[20,50],[20,40]]},
                               getdata=lambda key: arrays[key])

    def test_fallback_uses_saved_weights_masks_and_refined_parameter_count(self):
        stats = _gsas_fit_statistics(self.histogram(), {'varyList':['scale']})
        self.assertEqual(stats['chi2'], 10.0)
        self.assertEqual(stats['GoF'], round(math.sqrt(10),3))
        self.assertEqual(stats['Rwp'], round(100*math.sqrt(10/425),2))
        scaled = _gsas_fit_statistics(self.histogram(.25), {'varyList':['scale']})
        self.assertEqual(scaled['GoF'], round(math.sqrt(2.5),3))
        self.assertEqual(scaled['Rwp'], stats['Rwp'])

    def test_native_gof_is_reported_without_preset_rescaling(self):
        histogram = self.histogram()
        histogram.data['data'][0]['R'] = 5.339632
        stats = _gsas_fit_statistics(histogram, {
            'varyList':['scale'],
            'Rvals':{'Nvars':1, 'Rwp':5.716801, 'GOF':.75117006068}})
        self.assertEqual(stats, {'Rwp':5.72, 'Rp':5.34, 'chi2':.564, 'GoF':.751})


if __name__ == '__main__':
    unittest.main()
