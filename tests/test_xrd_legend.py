import copy
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from matplotlib.figure import Figure

from modules.xrd.xrd_plots import make_xrd_plot, normalize_legend_location


def fitted_pattern():
    tt = np.linspace(20, 80, 301)
    calculated = 10 + 100 * np.exp(-((tt - 55) / 1.5) ** 2)
    return {'tt': tt.tolist(), 'y_obs': calculated.tolist(),
            'y_calc': calculated.tolist(), 'y_background': [10.] * len(tt),
            'residuals': [0.] * len(tt), 'phase_results': [],
            'statistics': {'Rwp': 1, 'Rp': 1, 'chi2': 1, 'GoF': 1}}


class LegendRenderingTests(unittest.TestCase):
    def test_normalization_rejects_unsupported_locations(self):
        self.assertEqual(normalize_legend_location(), 'best')
        self.assertEqual(normalize_legend_location(' Outside Right '), 'outside right')
        with self.assertRaises(ValueError):
            normalize_legend_location('somewhere')

    def test_positions_render_inside_canvas_and_outside_clears_axes(self):
        result = fitted_pattern()
        original = copy.deepcopy(result)
        for location in ('best', 'upper left', 'lower right', 'outside right'):
            for theme in ('light', 'dark'):
                with self.subTest(location=location, theme=theme):
                    def inspect(figure, *_args, **_kwargs):
                        figure.canvas.draw()
                        axes = figure.axes[0]
                        legend = axes.get_legend()
                        bounds = legend.get_window_extent()
                        self.assertGreaterEqual(bounds.x0, 0)
                        self.assertGreaterEqual(bounds.y0, 0)
                        self.assertLessEqual(bounds.x1, figure.bbox.x1)
                        self.assertLessEqual(bounds.y1, figure.bbox.y1)
                        if location == 'outside right':
                            self.assertGreater(bounds.x0, axes.get_window_extent().x1)
                        else:
                            self.assertEqual(legend._loc, legend.codes[location])
                    with patch.object(Figure, 'savefig', inspect):
                        make_xrd_plot(result, {'legend_location': location},
                                      'unused.png', theme=theme)
        self.assertEqual(result['tt'], original['tt'])
        self.assertEqual(result['y_calc'], original['y_calc'])
        self.assertEqual(result['statistics'], original['statistics'])


class LegendApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import app
        cls.server = app

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.client = self.server.app.test_client()
        self.server._xrd_plot_cache.clear()
        self.paths = {theme: str(self.directory / f'{theme}.png')
                      for theme in ('light', 'dark')}
        for path in self.paths.values():
            Path(path).write_bytes(b'initial image')

    def test_new_fit_forwards_location_and_returns_regeneration_token_and_project(self):
        pattern = fitted_pattern()
        output = {'result': pattern, 'plot_path': self.paths['light'],
                  'plot_paths': self.paths, 'statistics': pattern['statistics'],
                  'phase_results': [], 'zero_shift': 0, 'summary_path': 'fit.xlsx',
                  'project_path': 'fit.gpx'}
        with patch.object(self.server, 'UPLOAD_DIR', str(self.directory)), \
                patch.object(self.server.xrd_processor, 'run', return_value=output) as run:
            response = self.client.post('/api/process_xrd', data={
                'file': (io.BytesIO(b'20 1\n21 2\n'), 'scan.xy'),
                'phases': json.dumps([{'name': 'test', 'cod_id': 'test'}]),
                'legend_location': 'outside right',
                'output_dir': str(self.directory),
            })
        self.assertEqual(response.status_code, 200, response.json)
        self.assertEqual(run.call_args.kwargs['params']['legend_location'], 'outside right')
        self.assertEqual(response.json['project_path'], 'fit.gpx')
        self.assertEqual(response.json['legend_location'], 'outside right')
        self.assertIn(response.json['plot_token'], self.server._xrd_plot_cache)

    def test_regeneration_updates_both_exports_without_refining(self):
        context = {'result': fitted_pattern(), 'metadata': {'sample_id': 'kept'},
                   'plot_paths': self.paths, 'plot_theme': 'light'}
        original = copy.deepcopy(context['result'])
        token = self.server._store_characterization_context(
            self.server._xrd_plot_cache, self.server._xrd_plot_cache_lock, context)
        with patch('modules.xrd.xrd_plots.make_xrd_plot') as render, \
                patch.object(self.server.xrd_processor, 'run') as refine:
            response = self.client.post('/api/xrd/regenerate_plot', json={
                'plot_token': token, 'legend_location': 'lower left'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json['legend_location'], 'lower left')
        self.assertEqual(render.call_count, 2)
        self.assertEqual({call.kwargs['theme'] for call in render.call_args_list},
                         {'light', 'dark'})
        for call in render.call_args_list:
            self.assertEqual(call.args[1]['legend_location'], 'lower left')
            self.assertEqual(call.args[1]['sample_id'], 'kept')
        refine.assert_not_called()
        self.assertEqual(context['result'], original)

    def test_expired_context_and_invalid_location_do_not_write_exports(self):
        context = {'result': fitted_pattern(), 'metadata': {},
                   'plot_paths': self.paths, 'plot_theme': 'light'}
        token = self.server._store_characterization_context(
            self.server._xrd_plot_cache, self.server._xrd_plot_cache_lock, context)
        with patch('modules.xrd.xrd_plots.make_xrd_plot') as render:
            missing = self.client.post('/api/xrd/regenerate_plot', json={
                'plot_token': 'expired', 'legend_location': 'lower left'})
            invalid = self.client.post('/api/xrd/regenerate_plot', json={
                'plot_token': token, 'legend_location': 'invalid'})
        self.assertEqual(missing.status_code, 410)
        self.assertEqual(invalid.status_code, 400)
        render.assert_not_called()


if __name__ == '__main__':
    unittest.main()
