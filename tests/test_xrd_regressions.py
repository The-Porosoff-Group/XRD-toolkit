import unittest
from pathlib import Path

from modules.xrd.cif_cache import mp_normal_cache_key, normalize_mp_id
from modules.xrd.cod_api import get_stick_pattern
from modules.xrd.crystallography import (
    expand_sites_from_cif,
    filter_reflections_by_relative_intensity,
    generate_reflections,
    parse_cif,
)
from modules.xrd.gsasii_backend import (
    _covariance_diagnostics,
    _prepared_cif_reference,
    _run_refinement_steps,
)
from modules.xrd.mp_api import (
    _fixture_cif_for,
    _parse as parse_mp_entries,
    _structure_dict_to_cif,
)


ROOT = Path(__file__).resolve().parents[1]


class FakeProject:
    def __init__(self):
        self.controls = []
        self.steps = []

    def set_Controls(self, control, value):
        self.controls.append((control, value))

    def do_refinements(self, steps):
        self.steps.append(steps)


class XrdRegressionTests(unittest.TestCase):
    @staticmethod
    def _ceo2_cif():
        from pymatgen.core import Lattice, Structure

        primitive = Structure(
            Lattice.rhombohedral(3.866071, 60.0),
            ["Ce", "O", "O"],
            [[0, 0, 0], [0.75, 0.75, 0.75], [0.25, 0.25, 0.25]],
        )
        return _structure_dict_to_cif(
            primitive.as_dict(),
            "mp-20194",
            "CeO2",
            {"number": 225, "symbol": "Fm-3m", "crystal_system": "Cubic"},
        )

    def test_cycles_are_project_controls_not_refinement_keys(self):
        project = FakeProject()
        _run_refinement_steps(project, [{
            "set": {"Background": {"refine": True}},
            "cycles": 11,
        }])

        self.assertEqual(project.controls, [("cycles", 11)])
        self.assertEqual(len(project.steps), 1)
        self.assertNotIn("cycles", project.steps[0][0])

    def test_covariance_reports_convergence_and_high_correlation(self):
        result = _covariance_diagnostics({
            "Rvals": {
                "converged": True,
                "DelChi2": 0.001,
                "Max shft/sig": 0.2,
            },
            "varyList": ["0:0:Scale", "0:0:Size;i"],
            "sig": [2.0, 5.0],
            "covMatrix": [[4.0, 9.8], [9.8, 25.0]],
        })

        self.assertTrue(result["converged"])
        self.assertEqual(len(result["high_correlations"]), 1)
        self.assertAlmostEqual(
            result["high_correlations"][0]["correlation_pct"], 98.0)

    def test_mp_1552_normal_fixture_is_conventional_pbcn(self):
        text = _fixture_cif_for("mp-1552", purpose="normal_import")
        self.assertIsNotNone(text)
        parsed = parse_cif(text)

        self.assertEqual(parsed.get("spacegroup_number"), 60)
        self.assertAlmostEqual(parsed.get("a"), 4.72854310, places=6)
        self.assertAlmostEqual(parsed.get("b"), 6.05260249, places=6)
        self.assertAlmostEqual(parsed.get("c"), 5.20975702, places=6)
        self.assertEqual(len(parsed.get("sites") or []), 2)

    def test_unsafe_raw_fixture_cannot_replace_normal_mp_import(self):
        self.assertIsNone(
            _fixture_cif_for("mp-1894", purpose="normal_import"))
        self.assertIsNotNone(_fixture_cif_for("mp-1894"))

    def test_prepared_reference_uses_exact_cif_setting(self):
        text = (ROOT / "fixtures" / "mo2c_pbcn_mp_1552.cif").read_text(
            encoding="utf-8")
        reference = _prepared_cif_reference(
            text, fallback={"a": 99, "b": 98, "c": 97})

        self.assertAlmostEqual(reference["a"], 4.72854310, places=6)
        self.assertAlmostEqual(reference["b"], 6.05260249, places=6)
        self.assertAlmostEqual(reference["c"], 5.20975702, places=6)
        self.assertEqual(reference["spacegroup_number"], 60)

    def test_mp_cache_key_is_versioned_and_role_specific(self):
        self.assertEqual(
            mp_normal_cache_key("mp-1552"), "mp:mp-1552:normal:v3")

    def test_current_mp_ids_decode_to_numeric_aliases(self):
        self.assertEqual(normalize_mp_id("mp-aaaabdws"), "mp-20194")
        self.assertEqual(normalize_mp_id("mp-aaaaadag"), "mp-2034")
        self.assertEqual(normalize_mp_id("mp-aaacfyxf"), "mp-1018659")
        self.assertEqual(normalize_mp_id("mp-1552"), "mp-1552")
        self.assertEqual(
            mp_normal_cache_key("mp-aaaaachs"),
            "mp:mp-1552:normal:v3",
        )

    def test_encoded_mp_search_id_routes_to_numeric_fixture(self):
        results = parse_mp_entries([{
            "material_id": "mp-aaaaachs",
            "formula_pretty": "Mo2C",
            "symmetry": {
                "symbol": "Pbcn",
                "number": 60,
                "crystal_system": "Orthorhombic",
            },
            "energy_above_hull": 0.0,
            "structure": {"lattice": {
                "a": 4.7, "b": 5.2, "c": 6.0,
                "alpha": 90, "beta": 90, "gamma": 90,
            }},
        }])

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]["mp_id"], "mp-1552")
        self.assertEqual(results[0]["mp_api_id"], "mp-aaaaachs")
        self.assertIn("_cif_text", results[0])
        self.assertEqual(results[0]["spacegroup_number"], 60)

    def test_ceo2_mp_conversion_is_conventional_and_asymmetric(self):
        text = self._ceo2_cif()
        parsed = parse_cif(text)

        self.assertEqual(parsed["spacegroup_number"], 225)
        self.assertEqual(parsed["spacegroup"], "Fm-3m")
        self.assertAlmostEqual(parsed["a"], 5.46745, places=4)
        self.assertAlmostEqual(parsed["b"], parsed["a"], places=6)
        self.assertAlmostEqual(parsed["c"], parsed["a"], places=6)
        self.assertEqual(len(parsed["sites"]), 2)
        self.assertIn("_space_group_symop_operation_xyz", text)

    def test_ceo2_preview_has_real_fluorite_ticks_without_forbidden_ticks(self):
        text = self._ceo2_cif()
        phase = parse_cif(text)
        phase["cif_text"] = text
        sticks = get_stick_pattern(phase, 1.54056, 20.0, 90.0)
        hkls = [stick["hkl"] for stick in sticks]

        self.assertEqual(hkls, [
            "(111)", "(200)", "(220)", "(311)", "(222)",
            "(400)", "(331)", "(420)", "(422)",
        ])
        self.assertNotIn("(100)", hkls)
        self.assertNotIn("(110)", hkls)
        strongest = max(sticks, key=lambda stick: stick["rel_int"])
        self.assertEqual(strongest["hkl"], "(111)")

    def test_ceo2_asymmetric_preview_matches_expanded_fit_ticks(self):
        text = self._ceo2_cif()
        parsed = parse_cif(text)
        args = (
            parsed["a"], parsed["b"], parsed["c"],
            parsed["alpha"], parsed["beta"], parsed["gamma"],
            parsed["system"], parsed["spacegroup_number"],
            1.54056, 20.0, 90.0,
        )
        preview = generate_reflections(
            *args, hkl_max=12, sites=parsed["sites"], site_policy="auto")
        expanded_sites = expand_sites_from_cif(text)
        fitted = generate_reflections(
            *args, hkl_max=12, sites=expanded_sites,
            site_policy="expanded_full_cell_sites")

        preview_ticks = filter_reflections_by_relative_intensity(preview)
        fitted_ticks = filter_reflections_by_relative_intensity(fitted)
        preview_signature = [
            (reflection[2], round(reflection[0], 6), round(reflection[3], 6))
            for reflection in preview_ticks
        ]
        fitted_signature = [
            (reflection[2], round(reflection[0], 6), round(reflection[3], 6))
            for reflection in fitted_ticks
        ]
        self.assertEqual(preview_signature, fitted_signature)

    def test_manual_cif_upload_uses_backend_preview_on_xrd_page(self):
        paths = (
            ROOT / "templates" / "index.html",
            ROOT / "templates" / "xrd_toolkit" / "index.html",
        )
        for path in (candidate for candidate in paths if candidate.exists()):
            template = path.read_text(encoding="utf-8")
            upload_start = template.index("async function uploadCIF")
            upload_end = template.index("function removePhase", upload_start)
            upload = template[upload_start:upload_end]
            self.assertIn("fetch('/api/xrd/preview_cif'", upload)
            self.assertLess(
                upload.index("fetch('/api/xrd/preview_cif'"),
                upload.index("const getVal ="),
            )
            self.assertIn("preview_unavailable: true", upload)
            self.assertIn("if (phase.preview_unavailable) return [];", template)

    def test_cif_parser_accepts_case_insensitive_and_dotted_tags(self):
        parsed = parse_cif("""data_test
_CELL_LENGTH_A 5.4
_CELL_LENGTH_B 5.4
_CELL_LENGTH_C 5.4
_SPACE_GROUP.IT_NUMBER 225
_SPACE_GROUP.NAME_H-M 'Fm-3m'
loop_
 _atom_site.label
 _atom_site.type_symbol
 _atom_site.fract_x
 _atom_site.fract_y
 _atom_site.fract_z
 ce1 ce 0 0 0
""")
        self.assertEqual(parsed["spacegroup_number"], 225)
        self.assertEqual(parsed["spacegroup"], "Fm-3m")
        self.assertEqual(parsed["sites"], [("Ce", 0.0, 0.0, 0.0, 1.0)])

    def test_w_wc_mo_moc_preview_regressions(self):
        cases = {
            "W": ({
                "formula": "W", "a": 3.1652, "b": 3.1652, "c": 3.1652,
                "alpha": 90, "beta": 90, "gamma": 90,
                "system": "cubic", "spacegroup_number": 229,
                "sites": [("W", 0, 0, 0, 1)],
            }, None, ["(110)", "(200)", "(211)", "(220)"]),
            "WC": ({"spacegroup_number": 187, "system": "hexagonal"},
                   "wc_p-6m2_mp_1894.cif",
                   ["(001)", "(100)", "(101)", "(210)", "(002)",
                    "(211)", "(200)", "(102)", "(201)"]),
            "Mo": ({"spacegroup_number": 229, "system": "cubic"},
                   "mo_metal_bcc_im3m.cif",
                   ["(110)", "(200)", "(211)", "(220)"]),
            "MoC": ({"spacegroup_number": 187, "system": "hexagonal"},
                    "moc_p-6m2_mp_2305.cif",
                    ["(001)", "(100)", "(101)", "(210)", "(002)",
                     "(211)", "(200)", "(102)", "(201)"]),
            "W2C": ({"spacegroup_number": 60, "system": "orthorhombic"},
                    "w2c_pbcn_mp_2034.cif",
                    ["(021)", "(002)", "(200)", "(121)", "(102)",
                     "(221)", "(202)", "(040)", "(023)", "(321)",
                     "(302)", "(042)", "(240)", "(223)", "(142)",
                     "(104)", "(400)", "(242)", "(204)"]),
            "Mo2C": ({"spacegroup_number": 60, "system": "orthorhombic"},
                     "mo2c_pbcn_mp_1552.cif",
                     ["(021)", "(002)", "(200)", "(121)", "(102)",
                      "(221)", "(202)", "(040)", "(023)", "(321)",
                      "(302)", "(042)", "(240)", "(223)", "(142)",
                      "(104)", "(400)", "(242)", "(204)"]),
        }
        for name, (overrides, fixture_name, expected) in cases.items():
            with self.subTest(name=name):
                phase = dict(overrides)
                if fixture_name:
                    text = (ROOT / "fixtures" / fixture_name).read_text(
                        encoding="utf-8")
                    phase = {**parse_cif(text), **phase, "cif_text": text}
                sticks = get_stick_pattern(phase, 1.54056, 20.0, 90.0)
                self.assertEqual(
                    [stick["hkl"] for stick in sticks], expected)


if __name__ == "__main__":
    unittest.main()
