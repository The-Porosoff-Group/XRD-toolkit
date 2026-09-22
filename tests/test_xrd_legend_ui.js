// Run with node tests/test_xrd_legend_ui.js; no browser dependencies.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.join(__dirname, '..');
const elements = new Map();
function element(id) {
  if (!elements.has(id)) elements.set(id, {value:'best', textContent:'', disabled:false});
  return elements.get(id);
}
let request, rendered;
const context = vm.createContext({
  document:{getElementById:element},
  xrdLastResult:{plot_token:'cached-fit', plot_data:{tt:[1, 2]}, statistics:{Rwp:2}},
  renderXrdPublicationPlot(data) { rendered = data; },
  fetch:async (url, options) => {
    request = {url, body:JSON.parse(options.body)};
    return {ok:true, json:async () => ({legend_location:request.body.legend_location})};
  },
});
vm.runInContext(fs.readFileSync(path.join(root, 'templates/xrd_shared/legend_settings.js'), 'utf8'), context);

(async () => {
  const position = context.xrdLegendLayout('outside right', {}, 2, 0.45);
  assert.ok(position.x > 1);
  assert.equal(position.xanchor, 'left');
  assert.equal(context.xrdLegendLayout('lower left', {}, 2, 0.45).yanchor, 'bottom');
  const automatic = context.xrdLegendLayout('best', {
    tt:[0, 0.8, 0.9, 1], y_obs:[1, 95, 100, 95], y_calc:[1, 95, 100, 95],
  }, 1, 0.45);
  assert.equal(automatic.xanchor, 'left', 'Automatic placement avoids the upper-right peak');

  await context.updateXrdLegendSelection({value:'outside right'});
  assert.equal(request.url, '/api/xrd/regenerate_plot');
  assert.equal(request.body.plot_token, 'cached-fit');
  assert.equal(request.body.legend_location, 'outside right');
  assert.equal(element('xrd-result-legend-location').value, 'outside right');
  assert.equal(element('xrd-legend-location').value, 'outside right');
  assert.match(element('xrd-legend-status').textContent, /both PNG exports/);
  assert.equal(rendered.statistics.Rwp, 2);
  assert.equal(element('xrd-result-legend-location').disabled, false);

  context.fetch = async () => ({ok:false, json:async () => ({error:'Fit expired'})});
  await context.updateXrdLegendSelection({value:'lower left'});
  assert.match(element('xrd-legend-status').textContent, /PNG exports not updated: Fit expired/);
  assert.equal(element('xrd-legend-location').disabled, false);

  for (const file of ['templates/index.html', 'templates/xrd_toolkit/index.html']) {
    if (!fs.existsSync(path.join(root, file))) continue;
    const html = fs.readFileSync(path.join(root, file), 'utf8');
    assert.match(html, /id="xrd-result-legend-location"/);
    assert.match(html, /fd.append\('legend_location'/);
    assert.match(html, /xrd_shared\/legend_settings.js/);
    assert.match(html, /GSAS-II Project \(GPX\)/);
  }
  console.log('XRD legend UI: positions, automatic peak avoidance, export updates, and failures passed.');
})().catch(error => { console.error(error); process.exitCode = 1; });
