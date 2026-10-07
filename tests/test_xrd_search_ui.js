// Run with node tests/test_xrd_search_ui.js; no third-party dependencies.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.join(__dirname, '..');
const elements = new Map();
function element(id) {
  if (!elements.has(id)) elements.set(id, {value: '', textContent: '', disabled: false});
  return elements.get(id);
}
let nextResponse = {}, lastRequest;
const context = vm.createContext({
  window: {mpKeySet: true},
  document: {getElementById: element, querySelector: () => ({value: 'mp'})},
  fetch: async (url, options) => {
    lastRequest = {url, ...options};
    return {ok: true, json: async () => nextResponse};
  },
});
vm.runInContext(fs.readFileSync(path.join(root, 'templates/xrd_shared/mp_key_settings.js'), 'utf8'), context);

(async () => {
  context.updateSourceBadge();
  assert.equal(element('mp-key-remove').disabled, false);
  assert.equal(element('mp-api-key').value, '');

  element('mp-api-key').value = 'replacement-key';
  nextResponse = {mp_key_set: true, message: 'Saved'};
  await context.manageMpKey('save');
  assert.equal(lastRequest.method, 'PUT');
  assert.equal(JSON.parse(lastRequest.body).api_key, 'replacement-key');
  assert.equal(element('mp-api-key').value, '');
  assert.equal(element('mp-key-feedback').textContent, 'Saved');

  nextResponse = {valid: true, message: 'Valid'};
  await context.manageMpKey('test');
  assert.equal(lastRequest.url, '/api/xrd/validate_mp_key');
  assert.equal(lastRequest.body, '{}');
  assert.equal(context.window.mpKeySet, true);

  nextResponse = {mp_key_set: false, message: 'Removed'};
  await context.manageMpKey('remove');
  assert.equal(lastRequest.method, 'DELETE');
  assert.equal(lastRequest.body, undefined);
  assert.equal(context.window.mpKeySet, false);
  assert.equal(element('mp-key-remove').disabled, true);
  assert.match(element('mp-key-badge').textContent, /no key/);

  context.fetch = async () => { throw new Error('Disconnected'); };
  element('mp-api-key').value = 'unsaved-new-key';
  await context.manageMpKey('save');
  assert.equal(element('mp-key-feedback').textContent, 'Disconnected');
  assert.equal(element('mp-api-key').value, 'unsaved-new-key');
  assert.equal(context.window.mpKeySet, false);
  assert.equal(element('mp-key-save').disabled, false);
  context.fetch = async () => ({ok: false, status: 405});
  await context.manageMpKey('save');
  assert.match(element('mp-key-feedback').textContent, /Restart the toolkit/);

  for (const file of ['templates/index.html', 'templates/xrd_toolkit/index.html']) {
    if (!fs.existsSync(path.join(root, file))) continue;
    const html = fs.readFileSync(path.join(root, file), 'utf8');
    const source = html.slice(html.indexOf('function filterCandidates()'), html.indexOf('function renderCandidates()'));
    const filterContext = vm.createContext({
      document: {getElementById: element}, renderCandidates() {},
      xrdCandidates: [
        {formula: 'W2C', name: 'W2C', search_text: 'w2c tungsten carbon', description: 'Hexagonal phase'},
        {formula: 'Fe', search_text: 'fe iron', description: 'Cubic phase'},
      ], xrdFiltered: [],
    });
    vm.runInContext(source, filterContext);
    for (const query of ['tungsten', 'TUNGSTEN hexagonal', '  tungsten  carbon ', 'W2C']) {
      element('xrd-filter').value = query;
      filterContext.filterCandidates();
      assert.equal(filterContext.xrdFiltered.length, 1, `${file}: ${query}`);
      assert.equal(filterContext.xrdFiltered[0].formula, 'W2C');
    }
    element('xrd-filter').value = '';
    filterContext.filterCandidates();
    assert.equal(filterContext.xrdFiltered.length, 2);
  }
  console.log('XRD UI: key replacement, testing, removal, failures, and keyword filters passed.');
})().catch(error => { console.error(error); process.exitCode = 1; });
