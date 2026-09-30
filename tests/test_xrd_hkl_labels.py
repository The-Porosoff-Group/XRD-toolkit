import math
from pathlib import Path
import unittest
from unittest.mock import patch

from modules.xrd.crystallography import (
    _display_hkl,
    d_spacing,
    generate_reflections,
    parse_cif,
)
from modules.xrd.presentation import format_hkl


# The monoclinic MP-559140 MoO2 cell and asymmetric coordinates used in the
# sample study. Explicit expansion keeps these tests independent of pymatgen.
MOO2_CELL = (5.54083333, 4.86854305, 5.64137481, 90.0, 120.175262, 90.0)
MOO2_ASYMMETRIC_SITES = (
    ('Mo', 0.77010223, 0.99217693, 0.78676165, 1.0),
    ('O', 0.88821594, 0.28097169, 0.62228793, 1.0),
    ('O', 0.60936583, 0.80211239, 0.40827061, 1.0),
)


def moo2_reflections():
    sites = []
    for element, x, y, z, occupancy in MOO2_ASYMMETRIC_SITES:
        for position in ((x, y, z), (-x, -y, -z),
                         (-x, y + 0.5, -z + 0.5),
                         (x, -y + 0.5, z + 0.5)):
            sites.append((element, *(value % 1 for value in position), occupancy))
    return generate_reflections(
        *MOO2_CELL, 'monoclinic', 14, 1.540593, 10.0, 90.0,
        hkl_max=8, sites=sites, site_policy='expanded_full_cell_sites',
    )


class HklLabelTests(unittest.TestCase):
    def test_non_cubic_labels_preserve_relative_signs(self):
        for system in ('monoclinic', 'triclinic', 'hexagonal', 'trigonal',
                       'tetragonal', 'orthorhombic'):
            with self.subTest(system=system):
                self.assertEqual(_display_hkl(1, 2, -3, system), (1, 2, -3))
                self.assertEqual(_display_hkl(-1, -2, 3, system), (1, 2, -3))
                self.assertEqual(_display_hkl(0, -2, 3, system), (0, 2, -3))
                self.assertEqual(_display_hkl(0, 0, -3, system), (0, 0, 3))
        self.assertEqual(_display_hkl(0, 0, 0, 'triclinic'), (0, 0, 0))

    def test_cubic_display_convention_is_unchanged(self):
        self.assertEqual(_display_hkl(-1, 2, -3, 'cubic'), (3, 2, 1))
        self.assertEqual(_display_hkl(0, -2, 0, 'CUBIC'), (2, 0, 0))

    def test_moo2_signed_families_have_distinct_correct_positions_and_labels(self):
        refs = moo2_reflections()
        selected = []
        for hkl in ((1, 1, -1), (1, 1, 1)):
            d = d_spacing(*hkl, *MOO2_CELL, 'monoclinic')
            position = math.degrees(2 * math.asin(1.540593 / (2 * d)))
            ref = min(refs, key=lambda value: abs(value[0] - position))
            self.assertAlmostEqual(ref[0], position, places=10)
            self.assertEqual(ref[2], hkl)
            selected.append(ref)
        self.assertGreater(abs(selected[0][0] - selected[1][0]), 5.0)
        self.assertNotEqual(format_hkl(selected[0][2]), format_hkl(selected[1][2]))
        self.assertEqual(format_hkl(selected[0][2]), '(1 1 1\u0305)')
        self.assertEqual(format_hkl(selected[1][2]), '(1 1 1)')
        # Every label must still describe the actual d-spacing of its line.
        for _, d, hkl, _ in refs:
            self.assertAlmostEqual(d_spacing(*hkl, *MOO2_CELL, 'monoclinic'),
                                   d, places=10)

    def test_label_fix_does_not_change_moo2_positions_intensities_or_count(self):
        corrected = moo2_reflections()
        with patch('modules.xrd.crystallography._display_hkl',
                   side_effect=lambda h, k, l, system: (abs(h), abs(k), abs(l))):
            previous = moo2_reflections()
        # This covers multiplicity and structure-factor weights as well as
        # d-spacing, ordering and reflection count; only labels may change.
        self.assertEqual([(r[0], r[1], r[3]) for r in corrected],
                         [(r[0], r[1], r[3]) for r in previous])

    def test_hexagonal_fixture_labels_match_their_actual_metric(self):
        fixtures = Path(__file__).resolve().parents[1] / 'fixtures'
        for filename in ('wc_p-6m2_mp_1894.cif', 'moc_p-6m2_mp_2305.cif'):
            with self.subTest(fixture=filename):
                phase = parse_cif((fixtures / filename).read_text(encoding='utf-8'))
                cell = tuple(phase[key] for key in
                             ('a', 'b', 'c', 'alpha', 'beta', 'gamma'))
                refs = generate_reflections(
                    *cell, 'hexagonal', 187, 1.54056, 20.0, 90.0, hkl_max=8)
                for _, d, hkl, _ in refs:
                    self.assertAlmostEqual(
                        d_spacing(*hkl, *cell, 'hexagonal'), d, places=10)
                by_hkl = {ref[2]: ref for ref in refs}
                for l in (0, 1):
                    # 2-1l has h^2+hk+k^2=3; erasing the minus sign gives7.
                    ref = by_hkl[(2, -1, l)]
                    expected_d = 1 / math.sqrt(4 / cell[0]**2 + l**2 / cell[2]**2)
                    self.assertAlmostEqual(ref[1], expected_d, places=10)
                    unsigned_d = d_spacing(2, 1, l, *cell, 'hexagonal')
                    unsigned_angle = math.degrees(2 * math.asin(1.54056 / (2 * unsigned_d)))
                    self.assertGreater(unsigned_angle, 90.0)


if __name__ == '__main__':
    unittest.main()
