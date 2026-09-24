// Run with node tests/test_xrd_instrument_ui.js.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const elements = new Map();
function makeElement() {
  return {value:'', textContent:'', checked:false, disabled:false, style:{}, files:[], children:[],
    replaceChildren(...children) { this.children = children; },
    append(child) { this.children.push(child); if (child.id) elements.set(child.id, child); },
    setAttribute() {}};
}
function element(id) { if (!elements.has(id)) elements.set(id, makeElement()); return elements.get(id); }
class Data { constructor() { this.fields = new Map(); } append(k,v) { this.fields.set(k,v); } get(k) { return this.fields.get(k); } }
const base = {id:'generic_flat_plate', label:'Generic flat plate', geometry:'bragg_brentano'};
const saved = {id:'local_test', label:'My bench', geometry:'bragg_brentano', local:true};
let requests = [];
let profiles = [base];
const context = vm.createContext({
  document: {getElementById:element, createElement:makeElement, addEventListener() {}},
  FormData:Data, xrdLastResult:null, downloadFile() {},
  fetch:async (url, request) => {
    requests.push({url, request});
    if (request?.method === 'POST') { profiles = [base, saved]; return {ok:true, json:async () => ({instrument:saved.id})}; }
    return {ok:true, json:async () => ({default:base.id, instruments:profiles})};
  },
});
vm.runInContext(fs.readFileSync('templates/xrd_shared/instrument_settings.js', 'utf8'), context);
(async () => {
  await context.refreshXrdInstruments();
  assert.equal(element('xrd-instrument').value, base.id);
  assert.equal(element('xrd-upload-profile').style.display, 'none');
  element('xrd-instrument').value = 'upload';
  context.xrdInstrumentChanged();
  assert.equal(element('xrd-upload-profile').style.display, '');
  assert.throws(() => context.appendXrdInstrumentSettings(new Data()), /Choose a .instprm/);
  element('xrd-instprm-file').files = [{name:'mine.instprm'}];
  element('xrd-upload-geometry').value = 'bragg_brentano';
  element('xrd-spectrum').value = 'single';
  let fd = new Data(); context.appendXrdInstrumentSettings(fd);
  assert.equal(fd.get('instprm_file').name, 'mine.instprm');
  assert.equal(fd.get('spectrum'), 'single');
  element('xrd-calibration-mode').checked = true;
  context.xrdCalibrationChanged();
  fd = new Data(); context.appendXrdInstrumentSettings(fd);
  assert.equal(fd.get('instprm_file'), undefined, 'calibration must start fresh');
  element('xrd-calibration-mode').checked = false;
  element('xrd-local-instrument-name').value = 'My bench';
  await context.saveXrdLocalInstrument(false);
  assert.equal(element('xrd-instrument').value, saved.id);
  assert.match(element('xrd-instrument-status').textContent, /saved locally and selected/);
  assert.equal(requests.find(r => r.request?.method === 'POST').request.body.get('geometry'), 'bragg_brentano');
  // Selecting a saved profile cannot accidentally submit a stale upload.
  fd = new Data(); context.appendXrdInstrumentSettings(fd);
  assert.equal(fd.get('instprm_file'), undefined);
  context.renderXrdCalibrationResult({calibration_token:'token', calibration_validation:{passed:false, reasons:['invalid width']}});
  assert.equal(element('xrd-save-calibration').disabled, true);
  context.xrdLastResult = {calibration_token:'valid'};
  context.renderXrdCalibrationResult({calibration_token:'valid', calibration_name:'Bench fit', calibration_validation:{passed:true}});
  element('xrd-calibration-mode').checked = true;
  await context.saveXrdLocalInstrument(true);
  assert.equal(element('xrd-calibration-mode').checked, false);
  assert.equal(element('xrd-instrument').value, saved.id);
  context.fetch = async () => ({ok:false, json:async () => ({error:'Invalid profile'})});
  element('xrd-local-instrument-name').value = 'Bad upload';
  await context.saveXrdLocalInstrument(false);
  assert.equal(element('xrd-instrument-status').textContent, 'Invalid profile');
  assert.equal(element('xrd-save-upload').disabled, false);
  console.log('Instrument UI: upload, geometry, persistence, stale-file isolation, calibration save, and errors passed.');
})().catch(error => { console.error(error); process.exitCode = 1; });
