"""Regressions for chemical-name search and local API-key management."""
import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import yaml
from modules.xrd import mp_api
from modules.xrd.search_terms import candidate_search_text, ELEMENT_NAMES


class ChemicalNameTests(unittest.TestCase):
    def test_names_are_resolved_before_short_formula_heuristic(self):
        queries = {
            'tungsten': ['W'], 'TUNGSTEN': ['W'], 'tin': ['Sn'],
            'carbon': ['C'], 'iron': ['Fe'], 'cerium': ['Ce'],
            'aluminum': ['Al'], 'sulphur': ['S'],
            'tungsten carbide': ['W', 'C'], 'W carbide': ['W', 'C'],
            'tungsten-carbide': ['W', 'C'], 'tungsten(VI) oxide': ['W', 'O'],
            'W C': ['W', 'C'], 'w c': ['W', 'C'],
        }
        for query, elements in queries.items():
            with self.subTest(query=query), patch.object(
                    mp_api, 'search_by_elements', return_value=[]) as search, patch.object(
                    mp_api, 'search_by_formula') as formula:
                mp_api.search_by_name(query, 'test-key', max_results=100, strict=False)
                search.assert_called_once_with(elements, 'test-key', strict=False,
                                               max_results=100, sort_by='formula')
                formula.assert_not_called()

    def test_formulas_still_work_in_name_field(self):
        for query, expected in [('W2C', 'W2C'), ('w2c', 'W2C'), ('wc', 'WC'),
                                ('Mo2C', 'Mo2C'), ('CeO2', 'CeO2'),
                                ('CO', 'CO'), ('NO', 'NO')]:
            with self.subTest(query=query), patch.object(
                    mp_api, 'search_by_formula', return_value=[]) as search:
                mp_api.search_by_name(query, 'test-key')
                search.assert_called_once_with(expected, 'test-key', 50, 'formula')

    def test_unsupported_text_is_not_silently_discarded(self):
        for query in ('quartz', 'cubic tungsten', 'tungsten banana'):
            with self.subTest(query=query), patch.object(mp_api, '_get') as request:
                result = mp_api.search_by_name(query, 'test-key')
                self.assertIn('description', result['error'])
                request.assert_not_called()

    def test_tungsten_uses_correct_api_filters(self):
        for strict, expected in [(True, {'chemsys': 'W', '_limit': 50}),
                                 (False, {'elements': 'W', '_limit': 50})]:
            with patch.object(mp_api, '_get', return_value=[]) as request:
                mp_api.search_by_name('tungsten', 'test-key', strict=strict)
                request.assert_called_once_with(expected, 'test-key')

    def test_candidate_index_contains_full_names_and_descriptions(self):
        text = candidate_search_text({'formula': 'W2C', 'description': 'Hexagonal phase',
                                      'mp_id': 'mp-1008625'})
        for term in ('tungsten', 'carbon', 'hexagonal', 'mp-1008625'):
            self.assertIn(term, text)
        self.assertNotIn('cobalt', text)
        self.assertEqual(len(ELEMENT_NAMES), 118)


class MpKeySettingsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import app
        cls.server = app

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.config_path = Path(self.temporary.name) / 'config.yaml'
        self.original = {'materials_project': {'api_key': 'old-test-key', 'extra': True},
                         'cache': {'max_size_mb': 321}, 'custom': {'value': 'keep'}}
        self.config_path.write_text(yaml.safe_dump(self.original), encoding='utf-8')
        for attribute, value in [('CONFIG_PATH', str(self.config_path)),
                                 ('CONFIG', copy.deepcopy(self.original)),
                                 ('MP_API_KEY', 'old-test-key')]:
            patcher = patch.object(self.server, attribute, value)
            patcher.start()
            self.addCleanup(patcher.stop)
        self.client = self.server.app.test_client()

    def test_replace_and_remove_persist_and_apply_without_restart(self):
        response = self.client.put('/api/xrd/mp_key', json={'api_key': '  replacement-key  '})
        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.json['mp_key_set'])
        self.assertNotIn('replacement-key', response.get_data(as_text=True))
        saved = yaml.safe_load(self.config_path.read_text(encoding='utf-8'))
        self.assertEqual(saved['materials_project']['api_key'], 'replacement-key')
        self.assertEqual(saved['cache'], self.original['cache'])
        self.assertEqual(saved['custom'], self.original['custom'])
        self.assertTrue(saved['materials_project']['extra'])
        self.assertEqual(self.server.MP_API_KEY, 'replacement-key')
        self.assertEqual(self.server.load_config()['materials_project']['api_key'], 'replacement-key')
        with patch.object(self.server, 'mp_search_name', return_value=[]) as search:
            self.client.post('/api/xrd/search', json={'source': 'mp', 'name': 'tungsten'})
            self.assertEqual(search.call_args.args[1], 'replacement-key')
        response = self.client.delete('/api/xrd/mp_key')
        self.assertEqual(response.status_code, 200)
        self.assertFalse(response.json['mp_key_set'])
        self.assertEqual(self.server.MP_API_KEY, '')
        self.assertEqual(self.server.load_config()['materials_project']['api_key'], '')
        self.assertFalse(self.client.get('/api/status').json['mp_key_set'])
        response = self.client.post('/api/xrd/search', json={'source': 'mp', 'name': 'tungsten'})
        self.assertEqual(response.status_code, 400)
        self.assertIn('API key', response.json['error'])

    def test_can_save_when_config_does_not_exist(self):
        self.config_path.unlink()
        response = self.client.put('/api/xrd/mp_key', json={'api_key': 'new-test-key'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(self.server.load_config()['materials_project']['api_key'], 'new-test-key')

    def test_invalid_requests_do_not_change_key(self):
        original = self.config_path.read_bytes()
        for payload in ({}, [], {'api_key': None}, {'api_key': 123},
                        {'api_key': ''}, {'api_key': 'a key with spaces'}):
            with self.subTest(payload=payload):
                response = self.client.put('/api/xrd/mp_key', json=payload)
                self.assertEqual(response.status_code, 400)
                self.assertEqual(self.config_path.read_bytes(), original)
                self.assertEqual(self.server.MP_API_KEY, 'old-test-key')

    def test_failed_write_preserves_file_and_active_key(self):
        original = self.config_path.read_bytes()
        with patch.object(self.server.os, 'replace', side_effect=OSError('write failed')):
            response = self.client.put('/api/xrd/mp_key', json={'api_key': 'replacement-key'})
        self.assertEqual(response.status_code, 500)
        self.assertEqual(self.config_path.read_bytes(), original)
        self.assertEqual(self.server.MP_API_KEY, 'old-test-key')
        self.assertEqual(list(self.config_path.parent.glob('*.tmp')), [])

    def test_bad_yaml_is_preserved_and_not_echoed(self):
        original = 'materials_project: [old-test-key'
        self.config_path.write_text(original, encoding='utf-8')
        response = self.client.delete('/api/xrd/mp_key')
        self.assertEqual(response.status_code, 500)
        self.assertNotIn('old-test-key', response.get_data(as_text=True))
        self.assertEqual(self.config_path.read_text(encoding='utf-8'), original)

    def test_test_key_does_not_save_and_can_use_existing_key(self):
        with patch.object(self.server, 'mp_validate_key', return_value=(True, 'Valid')) as validate:
            self.client.post('/api/xrd/validate_mp_key', json={})
            validate.assert_called_with('old-test-key')
            self.client.post('/api/xrd/validate_mp_key', json={'api_key': 'replacement-key'})
            validate.assert_called_with('replacement-key')
        self.assertEqual(self.server.MP_API_KEY, 'old-test-key')
        self.assertEqual(yaml.safe_load(self.config_path.read_text()), self.original)

    def test_views_include_controls_without_exposing_saved_key(self):
        for url in ('/', '/xrd'):
            with self.subTest(url=url):
                response = self.client.get(url)
                self.assertEqual(response.status_code, 200)
                html = response.get_data(as_text=True)
                self.assertIn('id="mp-key-remove"', html)
                self.assertIn('async function manageMpKey', html)
                self.assertNotIn('old-test-key', html)

    def test_keyword_priority_and_loose_mode_reach_search(self):
        with patch.object(self.server, 'mp_search_name', return_value=[]) as search:
            response = self.client.post('/api/xrd/search', json={
                'source': 'mp', 'elements': ['Fe'], 'name': 'tungsten', 'strict': False})
        self.assertEqual(response.status_code, 200)
        search.assert_called_once_with('tungsten', 'old-test-key',
                                       max_results=100, sort_by='formula', strict=False)


if __name__ == '__main__':
    unittest.main()
