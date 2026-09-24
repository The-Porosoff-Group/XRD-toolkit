import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from openpyxl import load_workbook

from modules.xrd import _write_summary_xlsx, run
from modules.xrd.fit_export import fit_parameter_rows
from test_xrd_publication_exports import publication_result


class FitParameterExportTests(unittest.TestCase):
    def test_native_values_survive_excel_precision_and_cell_limits(self):
        result = publication_result()
        cif = 'data_test\n' + 'atom 0.12345678901234567\n' * 2200
        result['gsas_native_parameters'] = {
            'final_project': {'data': {
                'Phases': {'W': {'Cell': [False, 3.1652345678901234],
                                 'Atoms': [['W1', 'W', '', 0.0, 0.0, 0.0]],
                                 'ranId': 1234567890123456789}},
                'Covariance': {'varyList': ['0:0:Scale'], 'sig': [0.000012345678901234]},
            }},
            'input_files': {'cif_text': cif},
        }
        result['fit_settings'] = {'zero': 0, 'enabled': False, 'absent': None,
                                  'special/key~': '=NOT_A_FORMULA()',
                                  'array': np.array([1.2345678901234567])}
        with tempfile.TemporaryDirectory() as tmp:
            path = _write_summary_xlsx(result, {'sample_id': 'W'}, 'GSAS-II', tmp)
            workbook = load_workbook(path)
            self.assertEqual(workbook.sheetnames, ['Summary', 'Plot Data', 'Fit Parameters'])
            sheet = workbook['Fit Parameters']
            groups = {}
            for row in sheet.iter_rows(min_row=2):
                section, pointer, kind, part, count, value = [cell.value for cell in row]
                self.assertEqual(row[-1].data_type, 's')
                self.assertLessEqual(len(value), 30000)
                groups.setdefault((section, pointer), []).append((part, value))
            decoded = {key: json.loads(''.join(value for _, value in sorted(parts)))
                       for key, parts in groups.items()}
            native = 'GSAS-II native parameters'
            self.assertEqual(decoded[(native, '/input_files/cif_text')], cif)
            self.assertEqual(decoded[(native, '/final_project/data/Phases/W/ranId')],
                             1234567890123456789)
            self.assertEqual(decoded[(native, '/final_project/data/Phases/W/Cell')],
                             [False, 3.1652345678901234])
            self.assertIs(decoded[('Fit settings', '/enabled')], False)
            self.assertIsNone(decoded[('Fit settings', '/absent')])
            self.assertEqual(decoded[('Fit settings', '/special~1key~0')], '=NOT_A_FORMULA()')
            workbook.close()

    def test_legacy_and_other_backends_do_not_claim_native_capture(self):
        legacy = fit_parameter_rows({}, {}, 'GSAS-II')
        other = fit_parameter_rows({}, {}, 'Le Bail')
        self.assertIn('Run the fit again', str(legacy))
        self.assertIn('Not applicable', str(other))

    def test_run_records_input_before_fit_and_exports_native_project(self):
        import base64
        result = publication_result()
        result['_gsas_project_base64'] = base64.b64encode(b'native-project-test').decode()
        phase = {'name': 'W', 'cif_text': 'data_W', 'a': 3.16}
        params = {'method': 'gsas2', 'phases': [phase], 'wavelength': 1.54056,
                  'instprm_file': 'test.instprm', 'legend_location': 'outside right',
                  'tt_min': 20, 'tt_max': 60, 'n_bg_coeffs': 'auto'}
        with tempfile.TemporaryDirectory() as tmp:
            scan = Path(tmp) / 'scan.xye'
            scan.write_text('20 100 2\n30 200 3\n60 100 2\n')
            profile = Path(tmp) / 'test.instprm'
            profile.write_text('Type:PXC\nLam:1.54056\nZero:0\nU:2\nV:-2\n'
                               'W:5\nX:0\nY:0\nSH/L:0.002\n')
            params['instprm_file'] = str(profile)
            with patch('modules.xrd.validate_phases', return_value=[phase]), \
                 patch('modules.xrd.gsasii_backend.is_available', return_value=True), \
                 patch('modules.xrd.gsasii_backend.run_gsas2', return_value=result), \
                 patch('modules.xrd.xrd_plots.make_xrd_plot') as plot:
                exported = run(str(scan), tmp, {'sample_id': 'W'}, params)
            settings = exported['result']['fit_settings']
            self.assertEqual(settings['parsed_input']['sigma'], [2, 3, 2])
            self.assertEqual(settings['effective_parameters']['iteration_limit'], 30)
            self.assertTrue(settings['effective_parameters']['auto_background'])
            phase['a'] = 99
            self.assertEqual(settings['input_phases'][0]['a'], 3.16)
            self.assertEqual(Path(exported['project_path']).read_bytes(), b'native-project-test')
            self.assertNotIn('_gsas_project_base64', exported['result'])
            self.assertEqual(plot.call_args.args[1]['legend_location'], 'outside right')


if __name__ == '__main__':
    unittest.main()
