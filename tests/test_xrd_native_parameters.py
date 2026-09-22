"""Native fit records retain reproducible values before display changes."""

import base64
import copy
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from modules.xrd.native_parameters import (
    NativeFitRecorder, native_value, project_snapshot, software_provenance,
)


class NativeVariable:
    def __init__(self):
        self.phase = 7
        self.name = 'Scale'

    def __str__(self):
        return '7:0:Scale'


def sample_project():
    return SimpleNamespace(
        names=[['Controls'], ['Phases', 'sample'], ['PWDR sample']],
        data={
            'Controls': {'data': {'max cyc': 20, 'deriv type': 'analytic Hessian'}},
            'Constraints': {'Phase': [[NativeVariable(), 1.0]]},
            'Covariance': {'data': {
                'varyList': ['0:0:Scale', ':0:Y'],
                'variables': np.array([0.12345678912345678, 0.8]),
                'sig': np.array([0.005, 0.02]),
                'covMatrix': np.array([[0.000025, 0.00001], [0.00001, 0.0004]]),
                'parmDict': {'0:0:Scale': 0.12345678912345678, ':0:SH/L': 0.002},
            }},
            'Phases': {'sample': {
                'General': {'Cell': [True, 3.1, 3.1, 3.1, 90., 90., 90., 29.791]},
                'Atoms': [['Ce1', 'Ce', 'U', 0., 0., 0., 1., 'm-3m', 4, 'I', 0.01]],
                'Histograms': {'PWDR sample': {
                    'Scale': [0.12345678912345678, True],
                    'Size': ['isotropic', [0.008, 1., 1.], [True, False, False]],
                    'Mustrain': ['isotropic', [1000., 1000., 1.], [False] * 3],
                    'Pref.Ori.': ['MD', 0.9, True, [0, 0, 1], 0, {}, []],
                    'new_future_parameter': [3.7, True],
                }},
            }},
            'PWDR sample': {
                'data': [{'wtFactor': 0.04}, [np.arange(5.) for _ in range(6)]],
                'Sample Parameters': {'Scale': [1., False], 'Shift': [2.5, True]},
                'Instrument Parameters': [{'U': [1.6, 1.6, False],
                                           'Y': [0., 0.8, True]}, {}],
                'Background': [['chebyschev-1', True, 2, 100., 0.2], {}],
                'Reflection Lists': {'sample': {'RefList': np.ones((4, 12))}},
            },
        },
    )


class NativeParameterCaptureTests(unittest.TestCase):
    def test_installed_gsas_version_identifies_source_commit_and_release(self):
        module = SimpleNamespace(GSASIIpath=SimpleNamespace(
            getSavedVersionInfo=lambda: SimpleNamespace(
                git_version='abc123', git_versiontag='v5.7.4', git_tags=['5833'])))
        software = software_provenance(module, __file__)
        self.assertEqual(software['gsas_version']['git_version'], 'abc123')
        self.assertEqual(software['gsas_version']['git_versiontag'], 'v5.7.4')
        self.assertEqual(software['gsas_version']['git_tags'], ['5833'])

    def test_all_native_fields_are_detached_without_a_parameter_allowlist(self):
        project = sample_project()
        notes = []
        snapshot = project_snapshot(project, notes)
        native = snapshot['data']
        self.assertEqual(native['PWDR sample']['Instrument Parameters'][0]['Y'],
                         [0., 0.8, True])
        hap = native['Phases']['sample']['Histograms']['PWDR sample']
        self.assertEqual(hap['new_future_parameter'], [3.7, True])
        self.assertEqual(hap['Size'][1][0], 0.008)
        self.assertEqual(native['Phases']['sample']['Atoms'][0][2], 'U')
        self.assertEqual(native['Covariance']['data']['varyList'],
                         ['0:0:Scale', ':0:Y'])
        self.assertEqual(native['Covariance']['data']['variables'][0],
                         0.12345678912345678)
        self.assertEqual(native['Covariance']['data']['covMatrix'][0][1], 0.00001)
        self.assertTrue(native['Phases']['sample']['General']['Cell'][0])
        project.data['Phases']['sample']['General']['Cell'][0] = False
        project.data['Covariance']['data']['sig'][0] = 1000
        self.assertTrue(native['Phases']['sample']['General']['Cell'][0])
        self.assertEqual(native['Covariance']['data']['sig'][0], 0.005)
        self.assertEqual(snapshot['excluded_paths'], [])
        self.assertTrue(notes)
        json.dumps(snapshot, allow_nan=False)

    def test_strict_json_preserves_masks_nonfinite_values_and_key_collisions(self):
        notes = []
        value = native_value({
            'nonfinite': [np.nan, np.inf, -np.inf],
            'mapping': {1: 'integer', '1': 'string'},
            'masked': np.ma.array([1., 2.], mask=[False, True]),
            'scalar': np.int64(20),
            'flag': np.bool_(True),
            'variable': NativeVariable(),
        }, notes)
        json.dumps(value, allow_nan=False)
        self.assertEqual(value['nonfinite'][0]['value'], 'nan')
        self.assertEqual(len(value['mapping']['entries']), 2)
        self.assertEqual(value['masked']['mask'], [False, True])
        self.assertEqual(value['masked']['data'], [1., 2.])
        self.assertEqual(value['scalar'], 20)
        self.assertIs(value['flag'], True)
        self.assertEqual(value['variable']['text'], '7:0:Scale')
        self.assertEqual(value['variable']['attributes']['phase'], 7)
        self.assertTrue(any('NativeVariable' in note for note in notes))

    def test_repeated_stage_snapshots_omit_only_named_derived_arrays(self):
        project = sample_project()
        original = copy.deepcopy(project.data)
        compact = project_snapshot(project, compact=True)
        self.assertNotIn('Covariance', compact['data'])
        self.assertNotIn('data', compact['data']['PWDR sample'])
        self.assertNotIn('Reflection Lists', compact['data']['PWDR sample'])
        self.assertEqual(len(compact['excluded_paths']), 3)
        self.assertEqual(compact['data']['PWDR sample']['Sample Parameters'],
                         original['PWDR sample']['Sample Parameters'])
        self.assertEqual(compact['data']['Phases']['sample']['Histograms'],
                         original['Phases']['sample']['Histograms'])
        self.assertIn('Covariance', project.data)
        self.assertIn('Reflection Lists', project.data['PWDR sample'])

    def test_final_snapshot_and_gpx_precede_destructive_display_changes(self):
        project = sample_project()
        recorder = NativeFitRecorder({'options': {'refine_cell': True}},
                                     SimpleNamespace(), __file__)
        stage = recorder.begin_stage(project, 'final constrained fit', 6,
                                     [{'set': {'Cell': True}, 'cycles': 15}])
        stage['outcome'] = {'failed': False}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'fit.gpx'
            project.save = lambda: path.write_bytes(b'native fitted project')
            recorder.finish(project, path)
            project.data['Phases']['sample']['General']['Cell'][0] = False
            project.data['Phases']['sample']['Histograms']['PWDR sample']['Scale'][0] = 0
            project.data['Controls']['data']['max cyc'] = 0
            project.data['Covariance']['data']['varyList'].clear()
            final = recorder.payload['final_project']['data']
            self.assertTrue(final['Phases']['sample']['General']['Cell'][0])
            self.assertEqual(final['Controls']['data']['max cyc'], 20)
            self.assertEqual(final['Phases']['sample']['Histograms']['PWDR sample']['Scale'][0],
                             0.12345678912345678)
            self.assertEqual(final['Covariance']['data']['varyList'],
                             ['0:0:Scale', ':0:Y'])
            self.assertEqual(stage['recipe'][0]['cycles'], 15)
            self.assertEqual(base64.b64decode(recorder.project_base64),
                             b'native fitted project')
        json.dumps(recorder.payload, allow_nan=False)

    def test_unavailable_native_archive_is_explicit_and_keeps_snapshot(self):
        recorder = NativeFitRecorder({}, SimpleNamespace(), __file__)
        recorder.finish(sample_project(), 'not-a-file.gpx')
        self.assertIsNone(recorder.project_base64)
        self.assertIn('final_project', recorder.payload)
        self.assertTrue(any('GPX companion could not be saved' in note
                            for note in recorder.payload['capture_notes']))


if __name__ == '__main__':
    unittest.main()
