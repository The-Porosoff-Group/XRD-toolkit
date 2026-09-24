"""Local profile persistence, geometry, and calibration route regressions."""
import io
import contextlib
import json
from pathlib import Path
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
from modules.xrd import instrument_profiles as profiles
from modules.xrd import gsasii_backend as backend
from modules.xrd import gsasii_calibration as calibration


PROFILE = '''#GSAS-II instrument parameter file
Type:PXC
Lam:1.540593
Zero:0
Polariz.:0.5
U:2
V:-2
W:5
X:0
Y:0
SH/L:0.002
Azimuth:0
'''


class InstrumentTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        patcher = patch.object(profiles, 'LOCAL_INSTRUMENT_DIR', self.root)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_unknown_scan_never_loads_smartlab(self):
        key, reason = backend.infer_instrument('Standard.txt')
        self.assertEqual(key, 'generic_flat_plate')
        profile = profiles.get_instrument_profile(key)
        self.assertEqual(profile['geometry'], 'bragg_brentano')
        self.assertIsNone(profiles.instrument_file(profile))
        self.assertNotIn('calibration_allow_x', profile)
        self.assertEqual(profiles.get_instrument_profile('generic'), profile)
        self.assertEqual(backend.infer_instrument('Rigaku_MiniFlex_Standard.txt')[0], key)
        self.assertEqual(backend.infer_instrument(metadata={'format': 'stepscan'})[0], key)
        with self.assertRaises(ValueError):
            profiles.get_instrument_profile('missing_local_profile')

    def test_save_reload_does_not_replace_previous_or_bundled_profile(self):
        original = dict(profiles.INSTRUMENT_PROFILES['smartlab'])
        first = profiles.save_local_instrument('Bench / Cu', 'bragg_brentano', PROFILE)
        second = profiles.save_local_instrument('Bench / Cu', 'capillary', PROFILE)
        self.assertNotEqual(first, second)
        saved = profiles.get_instrument_profiles()
        self.assertEqual(saved[first]['label'], 'Bench / Cu')
        self.assertEqual(saved[first]['displacement_param'], 'Shift')
        self.assertEqual(saved[second]['displacement_param'], 'DisplaceY')
        self.assertEqual(saved[first]['wavelength'], 1.540593)
        self.assertEqual(backend.get_instrument_profile(first), saved[first])
        self.assertEqual(Path(profiles.instrument_file(saved[first])).read_text(), PROFILE)
        self.assertEqual(profiles.INSTRUMENT_PROFILES['smartlab'], original)

    def test_primary_choices_exclude_legacy_geometry_and_saved_trials(self):
        local_key = profiles.save_local_instrument('Bench trial', 'bragg_brentano', PROFILE)
        primary = profiles.get_primary_instrument_profiles()
        self.assertEqual(list(primary), ['smartlab', 'synergy_s', 'benchtop_cu', 'none'])
        self.assertNotIn(local_key, primary)
        self.assertIn(local_key, profiles.get_instrument_profiles())
        for key in ('smartlab', 'synergy_s', 'benchtop_cu'):
            self.assertIsNotNone(profiles.instrument_file(primary[key]))

    def test_none_calibration_cannot_load_a_named_instrument_file(self):
        profile = profiles.get_instrument_profile('none')
        self.assertIsNone(profiles.instrument_file(profile))
        self.assertEqual(profile['geometry'], 'bragg_brentano')
        self.assertEqual(profile['zero_seed'], 0)
        self.assertTrue(profile['calibration'])
        self.assertNotIn('calibration_allow_x', profile)
        self.assertNotIn('calibration_allow_y', profile)
        self.assertNotIn('calibration_fixed_sh_l', profile)

    def test_reject_invalid_files_before_writing(self):
        invalid = ['', '# extra comment\n' + PROFILE, PROFILE.replace('PXC', 'PNT'), PROFILE.replace('U:2', 'U:nan'),
                   PROFILE.replace('Lam:1.540593', 'Lam:0'), PROFILE + 'U:1\n',
                   PROFILE.replace('SH/L:0.002', ''), PROFILE.replace('X:0', 'X:-1'),
                   PROFILE.replace('Polariz.:0.5', 'Polariz.:2')]
        for content in invalid:
            with self.subTest(content=content), self.assertRaises(ValueError):
                profiles.save_local_instrument('Bench', 'bragg_brentano', content)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_bundled_benchtop_profile_is_ready_without_local_upload(self):
        profile = profiles.get_instrument_profile('benchtop_cu')
        self.assertEqual(profile['geometry'], 'bragg_brentano')
        self.assertEqual(profile['calibration_range'], [20.0, 90.0])
        values = profiles.parse_instprm(Path(profiles.instrument_file(profile)).read_bytes())
        profiles.validate_profile_range(values, 20, 90)
        self.assertAlmostEqual(values['Lam1'], 1.540593)
        self.assertEqual(values['I(L2)/I(L1)'], 0.5)
        self.assertGreater(values['X'], 0)
        self.assertGreater(values['Y'], 0)

    def test_failed_metadata_write_leaves_no_half_saved_profile(self):
        with patch.object(profiles.os, 'replace', side_effect=OSError('disk full')):
            with self.assertRaises(OSError):
                profiles.save_local_instrument('Bench', 'bragg_brentano', PROFILE)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_geometry_is_explicit_for_single_wavelength_and_doublet(self):
        for geometry, expected in [('bragg_brentano', 'Bragg-Brentano'), ('capillary', 'Debye-Scherrer')]:
            histogram = types.SimpleNamespace(data={'Sample Parameters': {
                'Type': 'wrong', 'DisplaceY': [100, True], 'Shift': [20, True]}})
            profiles.configure_histogram_geometry(histogram, geometry)
            self.assertEqual(histogram.data['Sample Parameters']['Type'], expected)
            for key in ('Shift', 'DisplaceX', 'DisplaceY'):
                self.assertEqual(histogram.data['Sample Parameters'][key], [0, False])
            for key in ('Transparency', 'SurfRoughA', 'SurfRoughB', 'Absorption'):
                self.assertEqual(histogram.data['Sample Parameters'][key], [0, False])

    def test_fwhm_is_not_gaussian_sigma_and_full_range_is_checked(self):
        parameters = dict(U=0, V=0, W=100, X=0, Y=0)
        self.assertAlmostEqual(float(calibration._profile_fwhm_deg(parameters, [30])[0]),
                               np.sqrt(8 * np.log(2)) * 0.1)
        parameters.update(U=-6.4373, W=18.9295)
        self.assertTrue(calibration._profile_plausible(parameters, [20, 90])[0])
        self.assertFalse(calibration._profile_plausible(parameters, np.linspace(20, 150, 400))[0])
        parameters['U'] = float('nan')
        self.assertFalse(calibration._profile_plausible(parameters)[0])

    def test_single_wavelength_file_stays_single(self):
        with tempfile.TemporaryDirectory() as directory:
            path = backend._write_instprm(directory, 1.540593, kalpha2=False)
            data = profiles.parse_instprm(Path(path).read_bytes())
            self.assertEqual(data['Lam'], 1.540593)
            self.assertNotIn('Lam2', data)
            path = backend._write_instprm(directory, 1.540593, kalpha2=True)
            data = profiles.parse_instprm(Path(path).read_bytes())
            self.assertEqual(data['I(L2)/I(L1)'], 0.5)


class InstrumentRouteTests(InstrumentTests):
    @classmethod
    def setUpClass(cls):
        import app
        cls.server = app

    def setUp(self):
        super().setUp()
        self.client = self.server.app.test_client()
        patcher = patch.object(self.server, 'UPLOAD_DIR', str(self.root))
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_upload_roundtrip_and_both_interfaces(self):
        response = self.client.post('/api/xrd/instruments', data={
            'label': 'My bench', 'geometry': 'bragg_brentano',
            'instprm_file': (io.BytesIO(PROFILE.encode()), 'my.instprm')})
        self.assertEqual(response.status_code, 200, response.json)
        key = response.json['instrument']
        listed = self.client.get('/api/xrd/instruments').json
        self.assertEqual(listed['default'], 'generic_flat_plate')
        self.assertIn(key, [p['id'] for p in listed['instruments']])
        for route in ('/', '/xrd'):
            html = self.client.get(route).get_data(as_text=True)
            self.assertEqual(html.count('id="xrd-instrument"'), 1)
            self.assertNotIn('id="xrd-calibration-mode"', html)
            self.assertNotIn('id="xrd-local-instrument-name"', html)
            self.assertNotIn('id="xrd-spectrum"', html)
            self.assertIn('None / calibration', html)
            from toolkit_version import APP_VERSION
            self.assertIn(f'v{APP_VERSION}', html)

    def test_candidate_token_controls_geometry_and_range(self):
        token = self.server._store_characterization_context(
            self.server._xrd_calibration_cache, self.server._xrd_calibration_cache_lock,
            {'content': PROFILE, 'geometry': 'bragg_brentano',
             'fit_range': [20, 90], 'validation': {'passed': True}})
        response = self.client.post('/api/xrd/instruments', data={
            'label': 'Bench standard', 'calibration_token': token, 'geometry': 'capillary'})
        self.assertEqual(response.status_code, 200)
        profile = response.json['profile']
        self.assertEqual(profile['geometry'], 'bragg_brentano')
        self.assertEqual(profile['calibration_range'], [20, 90])
        self.assertEqual(self.client.post('/api/xrd/instruments', data={
            'label': 'Bad', 'calibration_token': 'unknown'}).status_code, 400)

    def test_calibration_uses_certified_si_without_database_phase(self):
        result = dict(
            tt=[28, 28.4, 29], y_obs=[1, 3, 1], y_calc=[1, 3, 1],
            y_background=[1, 1, 1], residuals=[0, 0, 0],
            statistics={'Rwp': 2}, params={'Zero': 0}, validation={'passed': True},
            geometry='bragg_brentano', fit_range=[20, 90], instprm_text=PROFILE,
            instprm_path=str(self.root / 'candidate.instprm'))
        def plot(result, metadata, path, **kwargs):
            Path(path).write_bytes(b'plot')
        with patch.object(calibration, 'run_calibration', return_value=result) as run, \
             patch('modules.xrd.xrd_plots.make_xrd_plot', side_effect=plot):
            response = self.client.post('/api/process_xrd', data={
                'file': (io.BytesIO(b'20 1\n21 2\n22 1\n'), 'standard.xy'),
                'method': 'gsas2', 'calibration_mode': 'true', 'phases': '[]',
                'instrument': 'generic_flat_plate', 'spectrum': 'single',
                'output_dir': str(self.root)})
        self.assertEqual(response.status_code, 200, response.json)
        kwargs = run.call_args.kwargs
        self.assertEqual(kwargs['phase']['formula'], 'Si')
        self.assertEqual(kwargs['phase']['a'], 5.431109)
        self.assertEqual(kwargs['instrument'], 'generic_flat_plate')
        self.assertEqual(kwargs['spectrum'], 'single')
        self.assertTrue(Path(kwargs['output_instprm']).is_relative_to(self.root))
        self.assertIn('calibration_token', response.json)

    def test_none_selection_alone_calibrates_with_no_uploaded_or_named_profile(self):
        result = dict(
            tt=[28, 28.4, 29], y_obs=[1, 3, 1], y_calc=[1, 3, 1],
            y_background=[1, 1, 1], residuals=[0, 0, 0],
            statistics={'Rwp': 2}, params={'Zero': 0}, validation={'passed': True},
            geometry='bragg_brentano', fit_range=[20, 90], instprm_text=PROFILE,
            instprm_path=str(self.root / 'candidate.instprm'))
        def plot(result, metadata, path, **kwargs):
            Path(path).write_bytes(b'plot')
        for geometry, instrument in [('bragg_brentano', 'none'), ('capillary', 'generic_capillary')]:
            with self.subTest(geometry=geometry), \
                 patch.object(calibration, 'run_calibration', return_value=dict(result, geometry=geometry)) as run, \
                 patch('modules.xrd.xrd_plots.make_xrd_plot', side_effect=plot):
                response = self.client.post('/api/process_xrd', data={
                    'file': (io.BytesIO(b'20 1\n21 2\n22 1\n'), 'standard.xy'),
                    'method': 'gsas2', 'phases': '[]', 'instrument': 'none',
                    'instrument_geometry': geometry,
                    'instprm_file': (io.BytesIO(b'invalid stale upload'), 'old.instprm'),
                    'output_dir': str(self.root)})
                self.assertEqual(response.status_code, 200, response.json)
                kwargs = run.call_args.kwargs
                self.assertEqual(kwargs['instrument'], instrument)
                profile = profiles.get_instrument_profile(kwargs['instrument'])
                self.assertIsNone(profiles.instrument_file(profile))
                self.assertEqual(profile['geometry'], geometry)
                self.assertEqual(kwargs['phase']['a'], 5.431109)
                self.assertEqual(kwargs['spectrum'], 'auto')
                self.assertEqual(list(self.root.glob('*.instprm')), [])
                self.assertIn('calibration_token', response.json)


@unittest.skipUnless(backend.is_available(), 'GSAS-II is not installed')
class CalibrationIntegrationTests(unittest.TestCase):
    def test_bundled_profile_imports_through_native_gsas_reader(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            instrument = profiles.instrument_file(profiles.get_instrument_profile('benchtop_cu'))
            expected = profiles.parse_instprm(Path(instrument).read_bytes())
            project = backend.G2sc.G2Project(newgpx=str(Path(directory) / 'bundled.gpx'))
            histogram = project.add_simulated_powder_histogram('standard', instrument, 20, 90, Tstep=0.02)
            actual = histogram.data['Instrument Parameters'][0]
            for key in ('Lam1', 'Lam2', 'I(L2)/I(L1)', 'U', 'V', 'W', 'X', 'Y', 'Zero', 'SH/L'):
                self.assertAlmostEqual(actual[key][1], expected[key])

    def test_single_wavelength_capillary_calibration_has_no_default_sample_strain(self):
        self._assert_synthetic_calibration('capillary', False)

    def test_single_wavelength_flat_plate_has_complete_sample_parameters(self):
        self._assert_synthetic_calibration('bragg_brentano', False)

    def test_doublet_capillary_has_complete_sample_parameters(self):
        self._assert_synthetic_calibration('capillary', True)

    def test_doublet_flat_plate_recovers_joint_lorentzian_broadening(self):
        self._assert_synthetic_calibration('bragg_brentano', True, mixed=True)

    def _assert_synthetic_calibration(self, geometry, doublet, mixed=False):
        # A known narrow profile would be absorbed by GSAS-II's default 1000
        # microstrain if the standard broadening were merely left unrefined.
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            root = Path(directory)
            phase = profiles.silicon_640g_phase()
            cif = root / 'standard.cif'
            cif.write_text(phase['cif_text'])
            inst = backend._write_instprm(str(root), 1.540593, kalpha2=doublet,
                                         u=0, v=0, w=4.6 if mixed else 10,
                                         x=2.1 if mixed else 0, y=6.9 if mixed else 0,
                                         sh_l=0.064 if mixed else 0.002)
            project = backend.G2sc.G2Project(newgpx=str(root / 'synthetic.gpx'))
            hist = project.add_simulated_powder_histogram('synthetic', inst, 20, 90, Tstep=0.02)
            profiles.configure_histogram_geometry(hist, geometry)
            standard = project.add_phase(str(cif), histograms=[hist])
            for hap in standard.data['Histograms'].values():
                hap['Scale'][0] = 1000.0
            standard.setSampleProfile(hist, 'size', 'isotropic', 10)
            standard.setSampleProfile(hist, 'microstrain', 'isotropic', 0)
            hist.data['Background'][0] = ['chebyschev-1', False, 1, 100]
            backend._run_refinement_steps(project, [{'set': {}, 'cycles': 0}])
            tt, y = hist.getdata('x'), hist.getdata('Ycalc')
            result = calibration.run_calibration(
                tt, y, np.sqrt(np.maximum(y, 1)), phase, 1.540593,
                tt_min=20, tt_max=90,
                instrument='generic_capillary' if geometry == 'capillary' else 'none',
                spectrum='cu_doublet' if doublet else 'single',
                output_instprm=str(root / 'calibrated.instprm'), keep_workdir=True)
            fitted = backend.G2sc.G2Project(result['project_path'])
            self.assertEqual(fitted.histograms()[0].data['Sample Parameters']['Type'],
                             'Debye-Scherrer' if geometry == 'capillary' else 'Bragg-Brentano')
            hap = next(iter(fitted.phases()[0].data['Histograms'].values()))
            self.assertEqual(hap['Mustrain'][1][0], 0)
            self.assertFalse(any(hap['Mustrain'][2]))
            self.assertEqual(hap['Size'][1][0], 10)
            self.assertFalse(any(hap['Size'][2]))
            self.assertFalse(fitted.phases()[0].data['General']['Cell'][0])
            values = profiles.parse_instprm(Path(result['candidate_instprm_path']).read_bytes())
            self.assertIn('Lam1' if doublet else 'Lam', values)
            self.assertEqual('Lam2' in values, doublet)
            self.assertLess(result['Rwp'], 5, json.dumps({
                'params': result['params'], 'stages': result['stage_log']}, indent=2))
            self.assertTrue(result['validation']['passed'])
            if mixed:
                self.assertEqual(result['lorentzian_term'], 'X+Y',
                                 json.dumps({'params': result['params'], 'Rwp': result['Rwp'],
                                             'stages': result['stage_log']}, indent=2))
                self.assertGreater(values['X'], 0)
                self.assertGreater(values['Y'], 0)


if __name__ == '__main__':
    unittest.main()
