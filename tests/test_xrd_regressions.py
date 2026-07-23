import unittest
from pathlib import Path

from modules.xrd.cif_cache import mp_normal_cache_key
from modules.xrd.crystallography import parse_cif
from modules.xrd.gsasii_backend import (
    _covariance_diagnostics,
    _prepared_cif_reference,
    _run_refinement_steps,
)
from modules.xrd.mp_api import _fixture_cif_for


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
            mp_normal_cache_key("mp-1552"), "mp:mp-1552:normal:v2")


if __name__ == "__main__":
    unittest.main()
