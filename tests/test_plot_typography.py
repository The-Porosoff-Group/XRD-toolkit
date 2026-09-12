"""Regressions for chemistry-aware labels, font consistency and text clipping."""

import unittest
from pathlib import Path
import subprocess
import sys

import matplotlib
from matplotlib.font_manager import FontProperties
from matplotlib.mathtext import MathTextParser
from PIL import Image

from modules.plot_style import (
    arial_plot, load_plot_font, plot_font_family, scientific_runs,
    scientific_mathtext, ScientificDraw, scientific_text_image,
)


class PlotTypographyTests(unittest.TestCase):
    def test_desktop_module_loading_can_import_shared_plot_style(self):
        # app.py discovers top-level modules, while tests import modules.*.
        # A clean interpreter catches imports that work in only one layout.
        root = Path(__file__).resolve().parents[1]
        result = subprocess.run([sys.executable, '-c',
            "import sys; sys.path.insert(0, 'modules'); "
            "import gc_processor; import xrd_processor"],
            cwd=root, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_formula_counts_do_not_change_sample_numbers_or_conditions(self):
        label = 'Sample 12: CO2 + 3 H2, 400 °C, 2026-09-12; W123'
        runs = scientific_runs(label)
        self.assertEqual([value for value, style in runs if style == 'sub'], ['2', '2'])
        self.assertEqual(''.join(value for value, _ in runs), label)
        self.assertEqual(scientific_runs('c2C4H8'),
                         [('c2C', 'normal'), ('4', 'sub'), ('H', 'normal'), ('8', 'sub')])
        self.assertIn(('0.8', 'sub'), scientific_runs('Ce0.8Zr0.2O2'))
        self.assertEqual(scientific_runs('WC1-x'), [('WC', 'normal'), ('1-x', 'sub')])
        self.assertEqual(scientific_mathtext(r'$\alpha$ at 400 °C'), r'$\alpha$ at 400 °C')

    def test_rotated_labels_include_scripts_and_descenders_inside_padding(self):
        for label in ('CO2 Conversion (%)', 'Quantity adsorbed (cm³ STP/g)',
                      '1 / [Q(p₀/p − 1)]', 'Fit (R²=0.99999)', 'CeO₂ (Fm3̅m)'):
            for size in (13, 28, 48):
                with self.subTest(label=label, size=size):
                    image = scientific_text_image(label, load_plot_font(size))
                    for rendered in (image, image.rotate(90, expand=True)):
                        left, top, right, bottom = rendered.getbbox()
                        self.assertGreaterEqual(left, 4)
                        self.assertGreaterEqual(top, 4)
                        self.assertLessEqual(right, rendered.width - 4)
                        self.assertLessEqual(bottom, rendered.height - 4)

    def test_measurement_contains_drawn_ink_for_centered_and_rotated_axis_anchors(self):
        for anchor in ('la', 'ma', 'ra', 'mm', 'ls', 'lt'):
            canvas = Image.new('RGBA', (900, 200))
            draw = ScientificDraw(canvas)
            font = load_plot_font(32)
            text = 'CO2 / H2, p/p₀, R²'
            measured = draw.textbbox((450, 100), text, font=font, anchor=anchor)
            draw.text((450, 100), text, font=font, fill='black', anchor=anchor)
            ink = canvas.getbbox()
            self.assertLessEqual(measured[0], ink[0])
            self.assertLessEqual(measured[1], ink[1])
            self.assertGreaterEqual(measured[2], ink[2])
            self.assertGreaterEqual(measured[3], ink[3])

    def test_math_and_pillow_use_same_font_and_export_restores_global_settings(self):
        previous = dict(matplotlib.rcParams)

        @arial_plot
        def check_font():
            formula = scientific_mathtext('CeO₂ (Fm3̅m)')
            parsed = MathTextParser('path').parse(
                formula + r' $R_{\mathrm{wp}}\ \chi^2$', prop=FontProperties(size=12))
            families = {font.family_name for font, *_ in parsed.glyphs}
            self.assertEqual(families, {plot_font_family()})
            self.assertEqual(load_plot_font(12).getname()[0], plot_font_family())

        check_font()
        self.assertEqual(dict(matplotlib.rcParams), previous)

        @arial_plot
        def failing_export():
            raise ValueError('deliberate export failure')

        with self.assertRaises(ValueError):
            failing_export()
        self.assertEqual(dict(matplotlib.rcParams), previous)


if __name__ == '__main__':
    unittest.main()
