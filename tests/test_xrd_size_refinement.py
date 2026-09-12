"""Regression coverage for broad peaks with the phase-card Size control."""
import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np

from modules.xrd import gsasii_backend as backend
from modules.xrd.crystallography import parse_cif


class SizeSeedTests(unittest.TestCase):
    def test_seed_tracks_nanocrystal_width_and_is_in_microns(self):
        tt = np.linspace(20, 90, 7001)
        wavelength = 1.54056
        expected_um = 0.008
        intensity = np.full_like(tt, 100.0)
        for center, height in ((28.65, 1000), (47.67, 500), (56.57, 400)):
            fwhm = np.degrees(wavelength / (
                10000 * expected_um * np.cos(np.radians(center / 2))))
            intensity += height * np.exp(-4 * np.log(2) * ((tt-center)/fwhm)**2)
        actual = backend._estimate_size_seed_um(tt, intensity, wavelength)
        self.assertAlmostEqual(actual, expected_um, delta=0.0001)

    def test_unusable_data_keeps_finite_default(self):
        for tt, y in (([1, 2], [1, 2]), (np.arange(10), np.ones(10)),
                      (np.arange(10), [np.nan] * 10)):
            with self.subTest(tt=tt):
                self.assertEqual(backend._estimate_size_seed_um(tt, y, 1.54), 1.0)


@unittest.skipUnless(backend.is_available(), 'GSAS-II is not installed')
class SizeRefinementIntegrationTests(unittest.TestCase):
    def _fit_synthetic(self, height_um=0.0, zero_mode=True,
                       geometry='bragg_brentano', start_a=5.46745):
        from pymatgen.core import Lattice, Structure
        from pymatgen.io.cif import CifWriter

        # Generate a known 8-nm pattern with GSAS-II itself. Start the fit
        # from a larger parent cell, as happens with a database CeO2 CIF.
        true_a = 5.391
        structure = Structure.from_spacegroup(
            'Fm-3m', Lattice.cubic(true_a), ['Ce', 'O'],
            [[0, 0, 0], [0.25, 0.25, 0.25]])
        true_cif = str(CifWriter(structure, symprec=0.01))
        structure.scale_lattice(start_a**3)
        start_cif = str(CifWriter(structure, symprec=0.01))
        phase = parse_cif(start_cif)
        phase.update(name='CeO2', system='cubic', Z=4, cif_text=start_cif)
        with tempfile.TemporaryDirectory() as tmp, \
                contextlib.redirect_stdout(io.StringIO()):
            inst = backend._write_instprm(
                tmp, 1.54056, u=1.6, v=-2.1, w=4.7, x=0, y=0,
                polariz=0.7, sh_l=0.002)
            cif_file = Path(tmp) / 'truth.cif'
            cif_file.write_text(true_cif, encoding='utf-8')
            project = backend.G2sc.G2Project(newgpx=str(Path(tmp) / 'truth.gpx'))
            histogram = project.add_simulated_powder_histogram(
                'synthetic', inst, 20, 90, Tstep=0.02)
            sample_type = ('Bragg-Brentano' if geometry == 'bragg_brentano'
                           else 'Debye-Scherrer')
            position_key = 'Shift' if geometry == 'bragg_brentano' else 'DisplaceY'
            histogram.data['Sample Parameters']['Type'] = sample_type
            histogram.data['Sample Parameters'][position_key] = [height_um, False]
            truth = project.add_phase(str(cif_file), histograms=[histogram])
            truth.set_HAP_refinements({'Size': {
                'type': 'isotropic', 'value': 0.008, 'refine': False}})
            histogram.data['Background'][0] = ['chebyschev-1', False, 1, 100.0]
            backend._run_refinement_steps(project, [{'set': {}, 'cycles': 0}])
            tt = histogram.getdata('x')
            y = histogram.getdata('Ycalc')
            result = backend.run_gsas2(
                tt, y, np.sqrt(np.maximum(y, 1)), [phase], 1.54056,
                tt_min=20, tt_max=90, n_bg_coeffs=6, max_cycles=20,
                auto_bg=False, instprm_file=inst, instrument='smartlab',
                options={
                    'geometry': geometry,
                    'verification_mode': True, 'verify_refine_cell': True,
                    'verify_use_zero_not_displace': zero_mode,
                    # A leftover Free X must not duplicate HAP Size.
                    'verify_refine_x': True, 'preferred_orientation': 'off',
                    'phase_options': [{'refine_cell': True, 'refine_size': True,
                                       'refine_mustrain': False, 'po_mode': 'off'}],
                })
        return result

    def test_size_recovers_broad_fluorite_peaks_with_offset_starting_cell(self):
        result = self._fit_synthetic()
        fitted = result['phase_results'][0]
        self.assertAlmostEqual(fitted['a'], 5.391, delta=0.002)
        self.assertAlmostEqual(fitted['crystallite_size_nm'], 8.0, delta=0.3)
        self.assertEqual(fitted['crystallite_size_source'], 'gsas_hap_size')
        self.assertEqual(fitted['gsas_size_model'], 'isotropic')
        self.assertEqual(fitted['gsas_size_lg_mix'], 1.0)
        self.assertAlmostEqual(fitted['gsas_hap_size_nm'], 8.0, delta=0.3)
        self.assertEqual(fitted['X'], 0.0)
        self.assertLess(abs(result['zero_shift']), 0.02)
        self.assertLess(result['statistics']['Rwp'], 2.0)
        self.assertFalse(result['refinement_diagnostics']['failed_stages'])
        self.assertEqual(result['displacement_param'], 'Shift')
        self.assertEqual(result['displacement_um'], 0.0)

    def test_flat_plate_height_recovers_known_offset_without_moving_zero(self):
        result = self._fit_synthetic(height_um=120.0, zero_mode=False)
        fitted = result['phase_results'][0]
        self.assertEqual(result['displacement_param'], 'Shift')
        self.assertAlmostEqual(result['displacement_um'], 120.0, delta=5.0)
        self.assertEqual(result['zero_shift'], 0.0)
        self.assertAlmostEqual(fitted['a'], 5.391, delta=0.0005)
        self.assertAlmostEqual(fitted['crystallite_size_nm'], 8.0, delta=0.3)
        self.assertLess(result['statistics']['Rwp'], 2.0)
        self.assertFalse(result['refinement_diagnostics']['failed_stages'])

    def test_capillary_displacement_uses_debye_scherrer_geometry(self):
        result = self._fit_synthetic(height_um=120.0, zero_mode=False,
                                     geometry='capillary', start_a=5.391)
        self.assertEqual(result['displacement_param'], 'DisplaceY')
        self.assertAlmostEqual(result['displacement_um'], 120.0, delta=5.0)
        self.assertEqual(result['zero_shift'], 0.0)
        self.assertAlmostEqual(result['phase_results'][0]['a'], 5.391, delta=0.0005)
        self.assertLess(result['statistics']['Rwp'], 2.0)


if __name__ == '__main__':
    unittest.main()
