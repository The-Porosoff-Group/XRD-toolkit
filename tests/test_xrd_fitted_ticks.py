"""Fitted tick display must not erase corrections or physical reflections."""
import copy
import unittest

from modules.xrd.gsasii_backend import _gsas_fitted_tick_refs


class FittedTickTests(unittest.TestCase):
    @staticmethod
    def reflection(hkl, pos, fc2, icorr):
        row = [0.0]*15
        row[:3] = hkl
        row[4:8] = [3.0, pos, 400.0, 100.0]
        row[9], row[11] = fc2, icorr
        return row

    def test_native_positions_and_signed_indices_survive_display_filter(self):
        rows = [self.reflection((1, 1, -1), 26.033, 20.0, 3.0),
                self.reflection((1, 1, 1), 36.765, 40.0, 2.0),
                self.reflection((0, 1, 0), 23.01, 0.01, 1.0)]
        record = {'RefList': rows, 'Super': False}
        before = copy.deepcopy(record)
        selected = _gsas_fitted_tick_refs(record, 10.0, 90.0)
        self.assertEqual([r[0] for r in selected], [26.033, 36.765])
        self.assertEqual([r[2] for r in selected], [(1, 1, -1), (1, 1, 1)])
        self.assertEqual(record, before)
        self.assertEqual(len(_gsas_fitted_tick_refs(
            record, 10.0, 90.0, include_weak=True)), 3)

    def test_filter_uses_predicted_intensity_and_never_modifies_native_list(self):
        # A weak structure factor need not imply a weak corrected intensity.
        rows = [self.reflection((1, 0, 0), 20., 100., .01),
                self.reflection((0, 1, 0), 30., .001, 1000.),
                self.reflection((0, 0, 1), 95., 100., 100.)]
        record = {'RefList': rows}
        selected = _gsas_fitted_tick_refs(record, 10., 90.)
        self.assertEqual(len(selected), 2)
        self.assertEqual(len(record['RefList']), 3)
        self.assertEqual(selected[0][3], selected[1][3])

    def test_superspace_ticks_preserve_signed_satellite_index(self):
        first = self.reflection((1, 1, -1), 26.033, 20.0, 3.0)
        second = self.reflection((1, 1, -1), 26.533, 20.0, 3.0)
        first.insert(3, -2)
        second.insert(3, 2)
        record = {'RefList': [first, second], 'Super': True}
        before = copy.deepcopy(record)
        selected = _gsas_fitted_tick_refs(record, 10., 90.)
        self.assertEqual([r[0] for r in selected], [26.033, 26.533])
        self.assertEqual([r[2] for r in selected],
                         [(1, 1, -1, -2), (1, 1, -1, 2)])
        self.assertEqual([r[3] for r in selected], [60., 60.])
        self.assertEqual(record, before)

    def test_malformed_display_records_warn_without_aborting_or_partial_ticks(self):
        valid = self.reflection((1, 1, 1), 36.7, 20., 3.)
        nonfinite = self.reflection((1, 1, 1), float('nan'), 20., 3.)
        invalid_index = self.reflection((1, 1, 1), 36.7, 20., 3.)
        invalid_index[0] = float('inf')
        for label, record in (
                ('missing record', None), ('missing list', {}),
                ('null list', {'RefList': None}),
                ('short row', {'RefList': [valid, [1., 2.]]}),
                ('nonfinite position', {'RefList': [valid, nonfinite]}),
                ('invalid index', {'RefList': [valid, invalid_index]})):
            with self.subTest(case=label), self.assertWarnsRegex(
                    UserWarning, 'ticks are unavailable.*fitted total is unchanged'):
                self.assertEqual(_gsas_fitted_tick_refs(record, 10., 90.), [])
        self.assertEqual(valid, self.reflection((1, 1, 1), 36.7, 20., 3.))


if __name__ == '__main__':
    unittest.main()
