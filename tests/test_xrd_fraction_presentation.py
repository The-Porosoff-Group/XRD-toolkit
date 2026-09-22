import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from modules.xrd import run
from modules.xrd.presentation import phase_legend_label
from scripts.xrd_batch import _compact_phase


ROOT = Path(__file__).resolve().parents[1]


class FractionPresentationTests(unittest.TestCase):
    def test_gsas_profile_seed_accepts_legacy_fits_without_mass_fractions(self):
        class BackendReached(Exception):
            pass

        phases = [{'name': 'phase', 'cif_text': 'provided CIF'}]
        seed = {
            'statistics': {'Rwp': 10},
            'phase_results': [
                {'weight_fraction_%': None, 'integrated_phase_fraction_%': 75,
                 'U': 0, 'V': 0, 'W': 0.001, 'X': 0, 'Y': 0.1},
                {'weight_fraction_%': None, 'integrated_phase_fraction_%': 25,
                 'U': 0, 'V': 0, 'W': 0.003, 'X': 0, 'Y': 0.3},
            ],
        }
        data = {'tt': np.array([20, 21]), 'intensity': np.array([1, 1]),
                'sigma': np.array([1, 1])}
        with tempfile.TemporaryDirectory() as output, \
                patch('modules.xrd.validate_phases', return_value=phases), \
                patch('modules.xrd.parse_xrd_file', return_value=data), \
                patch('modules.xrd.lebail.run_rietveld', return_value=seed), \
                patch('modules.xrd.gsasii_backend.is_available', return_value=True), \
                patch('modules.xrd.gsasii_backend.run_gsas2',
                      side_effect=BackendReached) as backend:
            source = Path(output) / 'pattern.xy'
            source.write_text('20 1\n21 1\n')
            with self.assertRaises(BackendReached):
                run(str(source), output, {}, {
                    'phases': phases, 'method': 'gsas2', 'instrument': 'generic'})
        self.assertAlmostEqual(backend.call_args.kwargs['seed_params']['W'], 15)
        self.assertAlmostEqual(backend.call_args.kwargs['seed_params']['Y'], 15)

    def test_batch_keeps_mass_and_diffraction_area_separate(self):
        phase = {
            'weight_fraction_%': 75.0,
            'weight_fraction_err_%': 1.2,
            'integrated_phase_fraction_%': 25.0,
            'weight_fraction_method': 'gsasii_mass',
            'weight_fraction_basis': 'modeled crystalline phases',
            'weight_fraction_note': 'Amorphous material excluded.',
        }
        row = _compact_phase('sample', {}, phase)
        self.assertEqual(row['weight_fraction_pct'], 75.0)
        self.assertEqual(row['weight_fraction_err_pct'], 1.2)
        self.assertEqual(row['diffraction_area_fraction_pct'], 25.0)
        self.assertEqual(row['integrated_phase_fraction_pct'], 25.0)
        self.assertEqual(row['weight_fraction_note'], phase['weight_fraction_note'])
        self.assertEqual(row['weight_fraction_basis'], phase['weight_fraction_basis'])
        phase['weight_fraction_%'] = None
        self.assertIsNone(_compact_phase('sample', {}, phase)['weight_fraction_pct'])

    def test_legend_never_substitutes_diffraction_area_for_mass(self):
        phase = {
            'formula': 'WC', 'weight_fraction_%': 75.0,
            'weight_fraction_err_%': 1.2, 'integrated_phase_fraction_%': 25.0,
        }
        self.assertEqual(phase_legend_label(phase), 'WC, 75.0 ± 1.2 wt. %')
        phase['weight_fraction_%'] = None
        self.assertEqual(phase_legend_label(phase), 'WC')

    @unittest.skipUnless(shutil.which('node'), 'Node.js required for UI checks')
    def test_both_interfaces_render_mass_and_preserve_missing_values(self):
        templates = [ROOT / 'templates/xrd_toolkit/index.html']
        if (ROOT / 'templates/index.html').exists():
            templates.append(ROOT / 'templates/index.html')
        for template in templates:
            with self.subTest(template=str(template)):
                source = template.read_text(encoding='utf-8')
                functions = []
                for name in ('fmtResultValue', '_fitNum', 'xrdWeightFractionHtml'):
                    match = re.search(
                        rf'^function {name}\(.*?^\}}', source, re.MULTILINE | re.DOTALL)
                    self.assertIsNotNone(match)
                    functions.append(match.group(0))
                script = '\n'.join(functions) + '''
const phase = {'weight_fraction_%': 75, 'weight_fraction_err_%': 1.2,
               'integrated_phase_fraction_%': 25};
console.log(JSON.stringify([
  xrdWeightFractionHtml(phase),
  xrdWeightFractionHtml({...phase, 'weight_fraction_%': null}),
  xrdWeightFractionHtml({...phase, 'weight_fraction_%': 0}),
  xrdWeightFractionHtml({...phase, 'weight_fraction_err_%': null}),
  _fitNum(null), _fitNum(undefined), _fitNum(''), _fitNum(0)
]));
'''
                completed = subprocess.run(
                    [shutil.which('node'), '-e', script], text=True,
                    capture_output=True, check=True, encoding='utf-8')
                self.assertEqual(json.loads(completed.stdout), [
                    '75 &plusmn; 1.2%', 'n/a', '0 &plusmn; 1.2%', '75%',
                    None, None, None, 0,
                ])
                self.assertIn('Diffraction area = </span>', source)
                self.assertIn('modeled crystalline phases', source)


if __name__ == '__main__':
    unittest.main()
