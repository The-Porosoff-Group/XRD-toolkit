// Run with node tests/test_plot_label_ui.js; exercises the actual preview formatter.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
for (const relative of ['../templates/index.html', '../templates/xrd_toolkit/index.html']) {
  const file = path.join(__dirname, relative);
  if (!fs.existsSync(file)) continue; // The standalone toolkit has only the XRD template.
  const html = fs.readFileSync(file, 'utf8');
  const context = vm.createContext({});
  for (const name of ['escHtml', 'scientificLabelHtml']) {
    const start = html.indexOf(`function ${name}(`);
    assert.ok(start >= 0, `${relative}: ${name}`);
    vm.runInContext(html.slice(start, html.indexOf('\n}', start) + 2), context);
  }
  assert.equal(context.scientificLabelHtml('Ce_{0.8}Zr_{0.2}O_{x}'),
    'Ce<sub>0.8</sub>Zr<sub>0.2</sub>O<sub>x</sub>');
  assert.equal(context.scientificLabelHtml('Fe^{3+} / m^{2} g^{-1}'),
    'Fe<sup>3+</sup> / m<sup>2</sup> g<sup>-1</sup>');
  assert.equal(context.scientificLabelHtml('sample_12 CeZrO_{x'), 'sample_12 CeZrO_{x');
  assert.equal(context.scientificLabelHtml('<img src=x onerror=alert(1)>_{<b>x</b>}'),
    '&lt;img src=x onerror=alert(1)&gt;<sub>&lt;b&gt;x&lt;/b&gt;</sub>');
}
console.log('Plot label UI: script codes render correctly and user text remains escaped.');
