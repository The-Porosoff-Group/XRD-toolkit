"""Capture a repeatable interface recipe without changing fit controls."""
import base64
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from modules.xrd import run
from test_xrd_publication_exports import publication_result
from toolkit_version import APP_VERSION


PROFILE = '''#GSAS-II instrument parameter file
Type:PXC
Lam:1.540593
Zero:0
U:2
V:-2
W:5
X:0
Y:0
SH/L:0.002
'''


class InputProvenanceTests(unittest.TestCase):
    def test_capture_preserves_original_files_and_controls_before_refinement(self):
        phase = {'name': 'W', 'cod_id': 'manual', 'source': 'manual',
                 'cif_filename': 'my tungsten.cif', 'cif_text': 'data_W', 'a': 3.16}
        controls = {'controls': {'wavelength': '1.54056', 'fix_y_value': '7.500',
                                'checkboxes': {'xrd-fix-y': False}},
                    'phase_options': [{'po_mode': 'off', 'po_axis': '0 0 0 1',
                                       'po_value': '0.905'}]}
        result = publication_result()
        result['_gsas_project_base64'] = base64.b64encode(b'native-project').decode()

        def refine(*args, **kwargs):
            phase['a'] = 99
            controls['controls']['fix_y_value'] = 'changed later'
            return result

        with tempfile.TemporaryDirectory() as tmp:
            scan = Path(tmp) / 'sanitized_scan.xy'
            scan.write_text('20 100 2\n30 200 3\n60 100 2\n')
            profile = Path(tmp) / 'random-upload.instprm'
            profile.write_text(PROFILE, encoding='utf-8')
            params = {'method': 'gsas2', 'phases': [phase], 'instrument': 'generic_flat_plate',
                      'wavelength': 1.54056, 'instprm_file': str(profile),
                      'instprm_original_filename': 'my bench.instprm',
                      'interface_settings': controls, 'n_bg_coeffs': 'auto',
                      'tt_min': 20, 'tt_max': 60}
            with patch('modules.xrd.validate_phases', return_value=[phase]), \
                    patch('modules.xrd.gsasii_backend.is_available', return_value=True), \
                    patch('modules.xrd.gsasii_backend.run_gsas2', side_effect=refine), \
                    patch('modules.xrd.xrd_plots.make_xrd_plot'), \
                    patch('modules.xrd._write_summary_xlsx', return_value='fit.xlsx'):
                exported = run(str(scan), tmp,
                               {'sample_id': 'W', 'source_file': 'Original scan.xy'}, params)
            settings = exported['result']['fit_settings']
            self.assertEqual(settings['software']['toolkit_version'], APP_VERSION)
            self.assertEqual(settings['source']['filename'], 'Original scan.xy')
            self.assertEqual(settings['source']['sha256'],
                             hashlib.sha256(scan.read_bytes()).hexdigest())
            self.assertEqual(settings['instrument_file'], {
                'filename': 'my bench.instprm',
                'sha256': hashlib.sha256(profile.read_bytes()).hexdigest()})
            self.assertEqual(settings['interface_settings']['controls']['fix_y_value'], '7.500')
            self.assertEqual(settings['interface_settings']['phase_options'][0]['po_axis'], '0 0 0 1')
            self.assertFalse(settings['interface_settings']['controls']['checkboxes']['xrd-fix-y'])
            self.assertEqual(settings['input_phases'][0]['a'], 3.16)
            self.assertEqual(settings['input_phases'][0]['cif_filename'], 'my tungsten.cif')
            self.assertEqual(settings['effective_parameters']['wavelength'], 1.540593)
            self.assertEqual(settings['project_file'], Path(exported['project_path']).name)
            self.assertEqual(Path(exported['project_path']).read_bytes(), b'native-project')

    def test_route_records_inactive_entries_without_enabling_them(self):
        import app as server
        snapshot = {'controls': {'instrument': 'upload', 'wavelength': '1.54056',
                                 'wavelength_source': 'Cu', 'fix_y_value': '7.500',
                                 'checkboxes': {'xrd-fix-y': False}},
                    'phase_options': [{'po_axis': '0 0 0 1', 'po_mode': 'off'}]}
        with tempfile.TemporaryDirectory() as tmp:
            plot = Path(tmp) / 'plot.png'
            plot.write_bytes(b'plot')
            output = {'result': {}, 'plot_path': str(plot),
                      'plot_paths': {'light': str(plot)}, 'statistics': {},
                      'phase_results': [], 'zero_shift': 0, 'summary_path': 'fit.xlsx'}
            with patch.object(server, 'UPLOAD_DIR', tmp), \
                    patch.object(server.xrd_processor, 'run', return_value=output) as refine:
                response = server.app.test_client().post('/api/process_xrd', data={
                    'file': (io.BytesIO(b'20 100\n30 200\n'), 'Original scan.xy'),
                    'phases': json.dumps([{'cod_id': 'manual', 'source': 'manual',
                                          'cif_filename': 'phase.cif', 'cif_text': 'data_W'}]),
                    'instprm_file': (io.BytesIO(PROFILE.encode()), 'my bench.instprm'),
                    'instrument': 'upload', 'method': 'gsas2', 'wavelength': '1.54056',
                    'interface_settings': json.dumps(snapshot), 'output_dir': tmp})
            self.assertEqual(response.status_code, 200, response.json)
            params = refine.call_args.kwargs['params']
            self.assertEqual(params['interface_settings']['controls'], snapshot['controls'])
            self.assertEqual(params['interface_settings']['phase_options'], snapshot['phase_options'])
            self.assertFalse(params['verify_fix_y'])
            self.assertIsNone(params['verify_y_fixed_value'])
            self.assertEqual(params['wavelength'], 1.540593)
            self.assertEqual(params['instprm_original_filename'], 'my bench.instprm')
            self.assertEqual(params['phases'][0]['cif_filename'], 'phase.cif')


if __name__ == '__main__':
    unittest.main()
