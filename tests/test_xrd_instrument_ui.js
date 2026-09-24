// Run with node tests/test_xrd_instrument_ui.js.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const elements = new Map();
function makeElement() {
  return {value:'', textContent:'', disabled:false, style:{}, files:[], children:[],
    replaceChildren(...children) { this.children = children; },
    append(child) { this.children.push(child); }, setAttribute() {}};
}
for (const id of ['xrd-instrument', 'xrd-instprm-file', 'xrd-instprm-override',
    'xrd-calibration-help', 'xrd-instrument-status', 'btn-xrd-gsas2',
    'xrd-calibration-phase-hint', 'xrd-calibration-geometry', 'xrd-calibration-result',
    'xrd-baseline-comparison', 'xrd-reusable-fit', 'xrd-phase-results', 'xrd-pymatgen-pill']) {
  elements.set(id, makeElement());
}
function element(id) { return elements.get(id) || null; }
class Data { constructor() { this.fields = new Map(); } append(k,v) { this.fields.set(k,v); } get(k) { return this.fields.get(k); } }
const context = vm.createContext({
  document: {getElementById:element, createElement:makeElement, addEventListener() {}},
  downloadFile() {},
});
vm.runInContext(fs.readFileSync('templates/xrd_shared/instrument_settings.js', 'utf8'), context);

context.restoreXrdInstrument({instrument:'benchtop_cu'});
assert.equal(context.isXrdCalibration(), false);
assert.equal(element('xrd-calibration-help').style.display, 'none');
assert.equal(element('btn-xrd-gsas2').textContent, 'GSAS-II Refinement');
let fd = new Data(); context.appendXrdInstrumentSettings(fd);
assert.equal(fd.get('instrument'), 'benchtop_cu');
assert.equal(fd.get('calibration_mode'), undefined);
assert.equal(fd.get('instprm_file'), undefined);

// Optional file override is confined to named sample instruments.
element('xrd-instprm-file').files = [{name:'mine.instprm'}];
fd = new Data(); context.appendXrdInstrumentSettings(fd);
assert.equal(fd.get('instprm_file').name, 'mine.instprm');
element('xrd-instrument').value = 'none';
context.xrdInstrumentChanged();
assert.equal(context.isXrdCalibration(), true);
assert.equal(element('xrd-instprm-file').disabled, true);
assert.equal(element('xrd-calibration-help').style.display, '');
assert.equal(element('btn-xrd-gsas2').textContent, 'Run calibration');
element('xrd-calibration-geometry').value = 'capillary';
fd = new Data(); context.appendXrdInstrumentSettings(fd);
assert.equal(fd.get('calibration_mode'), 'true');
assert.equal(fd.get('instrument_geometry'), 'capillary');
assert.equal(fd.get('calibration_standard'), 'Si640g');
assert.equal(fd.get('instprm_file'), undefined, 'calibration must ignore a stale file');
assert.equal(fd.get('spectrum'), undefined, 'source defaults replace the removed spectrum control');

context.restoreXrdInstrument({instrument:'smartlab', checkboxes:{'xrd-calibration-mode':true}});
assert.equal(element('xrd-instrument').value, 'none');
context.restoreXrdInstrument({instrument:'synergy_s', instrument_geometry:'bragg_brentano', checkboxes:{'xrd-calibration-mode':true}});
assert.equal(element('xrd-instrument').value, 'none');
assert.equal(element('xrd-calibration-geometry').value, 'capillary', 'legacy named calibration keeps its measurement geometry');
context.restoreXrdInstrument({instrument:'generic_capillary'});
assert.equal(element('xrd-instrument').value, 'none');
assert.equal(element('xrd-calibration-geometry').value, 'capillary');
assert.match(element('xrd-instrument-status').textContent, /older profile/);
context.restoreXrdInstrument({instrument:'synergy_s'});
assert.equal(element('xrd-instrument').value, 'synergy_s');
assert.equal(element('xrd-instprm-file').value, '');
assert.equal(element('xrd-instprm-file').disabled, false);
element('xrd-instrument').value = 'local_trial';
assert.throws(() => context.appendXrdInstrumentSettings(new Data()), /Select an instrument/);

context.renderXrdCalibrationResult({calibration_token:'token', calibration_validation:{passed:false, reasons:['invalid width']}, candidate_instprm_path:'candidate.instprm'});
assert.match(element('xrd-calibration-result').children[0].textContent, /checks failed/);
assert.ok(element('xrd-calibration-result').children.some(e => e.textContent === 'Download candidate .instprm'));
context.renderXrdCalibrationResult({});
assert.equal(element('xrd-calibration-result').style.display, 'none');

const controls = fs.readFileSync('templates/xrd_shared/instrument_controls.html', 'utf8');
assert.deepEqual([...controls.matchAll(/<option value="([^"]+)"/g)].map(match => match[1]),
  ['smartlab', 'synergy_s', 'benchtop_cu', 'none']);
for (const file of ['templates/index.html', 'templates/xrd_toolkit/index.html']) {
  if (!fs.existsSync(file)) continue;
  const html = fs.readFileSync(file, 'utf8');
  const selector = html.indexOf("include 'xrd_shared/instrument_controls.html'");
  assert.ok(selector > html.indexOf('id="btn-xrd-gsas2"') && selector < html.indexOf('class="xrd-preset-row"'));
  assert.doesNotMatch(html, /id="xrd-(calibration-mode|local-instrument-name|spectrum|legend-location)"/);
}
console.log('Instrument UI: four choices, calibration dispatch, file isolation, legacy presets, results, and restored placement passed.');
