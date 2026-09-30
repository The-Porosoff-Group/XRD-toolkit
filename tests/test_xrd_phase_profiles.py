"""Native component curves must preserve fitted broadening and intensities."""
import contextlib
import copy
import io
import pickle
from pathlib import Path
import tempfile
import unittest

import numpy as np

from modules.xrd import gsasii_backend as backend


@unittest.skipUnless(backend.is_available(), 'GSAS-II is not installed')
class NativePhaseProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pymatgen.core import Lattice, Structure
        from pymatgen.io.cif import CifWriter

        # Expected components come from GSAS-II's own PhasePartials output,
        # not a second implementation of the reconstruction helper. This
        # calculates fixed parameters with zero cycles; it does not fit data.
        carbide = Structure.from_spacegroup(
            'Fm-3m', Lattice.cubic(4.18), ['Mo', {'C': 0.5}],
            [[0, 0, 0], [0.5, 0.5, 0.5]])
        oxide = Structure.from_spacegroup(
            136, Lattice.tetragonal(4.8, 3.1), ['Mo', 'O'],
            [[0, 0, 0], [0.3, 0.3, 0]])
        with tempfile.TemporaryDirectory() as directory, \
                contextlib.redirect_stdout(io.StringIO()):
            folder = Path(directory)
            instrument = backend._write_instprm(
                directory, 1.540593, kalpha2=True, polariz=0.63,
                sh_l=0.03, u=8.0, v=0.0, w=20.0, x=0.0, y=0.0,
                zero_seed=-0.12)
            project = backend.G2sc.G2Project(
                newgpx=str(folder / 'components.gpx'))
            histogram = project.add_simulated_powder_histogram(
                'synthetic components', instrument, 10, 90, Tstep=0.02)
            histogram.data['Sample Parameters']['Type'] = 'Debye-Scherrer'
            histogram.data['Sample Parameters']['Scale'] = [3.0, False]
            histogram.data['Background'][0] = [
                'chebyschev-1', False, 1, 25.0]
            # Avoid Poisson-generated observations in a dummy calculation.
            histogram.data['data'][1][1][:] = 100.0
            inst = histogram.data['Instrument Parameters'][0]
            inst['I(L2)/I(L1)'][1] = 0.37
            # Initial values are intentionally different from current ones.
            # Native GSAS uses index 1 for the actual calculation.
            inst['Lam1'][0] = 0.71073
            inst['Lam2'][0] = 0.71359
            inst['I(L2)/I(L1)'][0] = 0.1
            inst['SH/L'][0] = 0.002
            for name, structure, size, mix, scale in (
                    ('broad carbide', carbide, 0.0016, 0.4, 2.0),
                    ('narrow oxide', oxide, 0.035, 0.85, 0.7)):
                cif = folder / (name.replace(' ', '_') + '.cif')
                cif.write_text(str(CifWriter(structure, symprec=0.001)),
                               encoding='utf-8')
                phase = project.add_phase(str(cif), phasename=name,
                                          histograms=[histogram])
                phase.set_HAP_refinements({'Size': {
                    'type': 'isotropic', 'value': size, 'refine': False,
                    'LGmix': {'value': mix, 'refine': False}}})
                phase.data['Histograms'][histogram.name]['Scale'] = [
                    scale, False]
            partials_file = folder / 'native_partials.pickle'
            project.data['Controls']['data']['PhasePartials'] = str(partials_file)
            project.data['Controls']['data']['deriv type'] = 'analytic Hessian'
            project.set_Controls('cycles', 0)
            project.refine()
            histogram = project.histograms()[0]
            cls.tt = np.asarray(histogram.getdata('x')).copy()
            cls.total = (np.asarray(histogram.getdata('ycalc')) -
                         np.asarray(histogram.getdata('background')))
            cls.reflections = copy.deepcopy(histogram.data['Reflection Lists'])
            cls.instrument = copy.deepcopy(
                histogram.data['Instrument Parameters'][0])
            cls.phase_names = [phase.name for phase in project.phases()]
            cls.native_components = {}
            with partials_file.open('rb') as stream:
                if pickle.load(stream) is not None:
                    raise AssertionError('Missing native histogram header')
                pickle.load(stream)  # histogram id
                native_x = pickle.load(stream)
                pickle.load(stream)  # background
                np.testing.assert_array_equal(native_x, cls.tt)
                for _ in cls.phase_names:
                    name = pickle.load(stream)
                    cls.native_components[name] = pickle.load(stream)

    def reconstruct(self, reflections=None, instrument=None, names=None):
        return backend._compute_gsas_cw_phase_profiles(
            self.tt,
            self.reflections if reflections is None else reflections,
            self.instrument if instrument is None else instrument,
            self.phase_names if names is None else names)

    def test_broad_and_narrow_components_match_native_calculation(self):
        actual = self.reconstruct()
        self.assertEqual(len(actual), 2)
        for name, profile in zip(self.phase_names, actual):
            with self.subTest(phase=name):
                self.assertIsInstance(profile, np.ndarray)
                np.testing.assert_allclose(
                    profile, self.native_components[name], rtol=2e-12,
                    atol=2e-10)
        np.testing.assert_allclose(np.sum(actual, axis=0), self.total,
                                   rtol=2e-12, atol=2e-10)
        # The two native phases actually exercise different peak widths.
        broad = self.reflections[self.phase_names[0]]['RefList']
        narrow = self.reflections[self.phase_names[1]]['RefList']
        self.assertGreater(np.median(broad[:, 6]),
                           10 * np.median(narrow[:, 6]))

    def test_doublet_and_current_instrument_values_are_preserved(self):
        reference = self.reconstruct()
        changed_initial = copy.deepcopy(self.instrument)
        for key in ('Lam1', 'Lam2', 'I(L2)/I(L1)', 'SH/L'):
            changed_initial[key][0] *= 1.7
        for a, b in zip(reference, self.reconstruct(instrument=changed_initial)):
            np.testing.assert_array_equal(a, b)
        single = copy.deepcopy(self.instrument)
        single['I(L2)/I(L1)'][1] = 0.0
        for both, first in zip(reference, self.reconstruct(instrument=single)):
            self.assertGreater(float(np.linalg.norm(both - first)), 0.0)

    def test_superspace_column_offset_and_requested_phase_order(self):
        shifted = copy.deepcopy(self.reflections)
        for entry in shifted.values():
            entry['RefList'] = np.insert(entry['RefList'], 3, 0.0, axis=1)
            entry['Super'] = True
        actual = self.reconstruct(reflections=shifted,
                                  names=list(reversed(self.phase_names)))
        for name, profile in zip(reversed(self.phase_names), actual):
            np.testing.assert_allclose(profile, self.native_components[name],
                                       rtol=2e-12, atol=2e-10)

    def test_reconstruction_does_not_mutate_native_inputs(self):
        original = pickle.dumps((self.tt, self.reflections, self.instrument),
                                protocol=pickle.HIGHEST_PROTOCOL)
        self.reconstruct()
        self.assertEqual(original, pickle.dumps(
            (self.tt, self.reflections, self.instrument),
            protocol=pickle.HIGHEST_PROTOCOL))

    def test_one_point_endpoint_window_preserves_its_only_sample(self):
        # GSAS's profile wrapper needs two samples for an auxiliary integral.
        # A clipped one-sample tail still has a valid pointwise intensity.
        grid = np.array([20.0, 21.0, 22.0])
        instrument = copy.deepcopy(self.instrument)
        instrument['SH/L'][1] = 0.002
        instrument['I(L2)/I(L1)'][1] = 0.0
        for index in (0, 2):
            with self.subTest(endpoint=index):
                row = np.zeros((1, 15))
                row[0, [5, 6, 7, 9, 11]] = [grid[index], 1e-8, 0, 2, 3]
                actual = backend._compute_gsas_cw_phase_profiles(
                    grid, {'edge': {'RefList': row}}, instrument, ['edge'])[0]
                native = backend.G2pwd.getFCJVoigt3(
                    grid[index], 1e-8, 0, 0.002, grid)[0]
                self.assertAlmostEqual(actual[index], 6 * native[index])
                np.testing.assert_array_equal(
                    np.delete(actual, index), np.zeros(2))

    def test_missing_or_malformed_reflections_are_rejected(self):
        missing_phase = copy.deepcopy(self.reflections)
        del missing_phase[self.phase_names[0]]
        missing_list = copy.deepcopy(self.reflections)
        del missing_list[self.phase_names[0]]['RefList']
        short_rows = copy.deepcopy(self.reflections)
        short_rows[self.phase_names[0]]['RefList'] = np.ones((1, 7))
        bad_width = copy.deepcopy(self.reflections)
        bad_width[self.phase_names[0]]['RefList'][0, 6] = -1.0
        nonfinite = copy.deepcopy(self.reflections)
        nonfinite[self.phase_names[0]]['RefList'][0, 9] = np.nan
        for label, reflections in (
                ('missing phase', missing_phase), ('missing list', missing_list),
                ('short rows', short_rows), ('negative variance', bad_width),
                ('nonfinite structure factor', nonfinite)):
            with self.subTest(case=label), self.assertRaises(ValueError):
                self.reconstruct(reflections=reflections)


if __name__ == '__main__':
    unittest.main()
