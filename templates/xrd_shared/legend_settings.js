function xrdLegendLayout(location, plotData, phaseCount, mainBottom) {
  const positions = {
    'upper right': {x:0.992, y:0.985, xanchor:'right', yanchor:'top'},
    'upper left': {x:0.008, y:0.985, xanchor:'left', yanchor:'top'},
    'lower right': {x:0.992, y:mainBottom + 0.015, xanchor:'right', yanchor:'bottom'},
    'lower left': {x:0.008, y:mainBottom + 0.015, xanchor:'left', yanchor:'bottom'},
    'center right': {x:0.992, y:(1 + mainBottom) / 2, xanchor:'right', yanchor:'middle'},
    'center left': {x:0.008, y:(1 + mainBottom) / 2, xanchor:'left', yanchor:'middle'},
    'upper center': {x:0.5, y:0.985, xanchor:'center', yanchor:'top'},
    'lower center': {x:0.5, y:mainBottom + 0.015, xanchor:'center', yanchor:'bottom'},
    'center': {x:0.5, y:(1 + mainBottom) / 2, xanchor:'center', yanchor:'middle'},
    'outside right': {x:1.02, y:0.985, xanchor:'left', yanchor:'top'},
  };
  if (positions[location]) return positions[location];
  // Plotly has no automatic placement. Score the data beneath each available
  // inside position; the PNG renderer uses Matplotlib's native best placement.
  const tt = plotData.tt || [];
  const observed = plotData.y_obs || [];
  const calculated = plotData.y_calc || [];
  let maxY = 1;
  for (const values of [observed, calculated]) {
    for (const value of values) maxY = Math.max(maxY, Number(value) || 0);
  }
  const span = tt.length > 1 ? tt[tt.length - 1] - tt[0] : 1;
  const width = 0.35;
  const height = Math.min(0.5, (phaseCount + 3) * 0.042);
  let best = positions['upper right'], bestScore = Infinity;
  for (const [name, position] of Object.entries(positions)) {
    if (name === 'outside right') continue;
    const left = position.x - (position.xanchor === 'right' ? width : position.xanchor === 'center' ? width / 2 : 0);
    const bottom = position.y - (position.yanchor === 'top' ? height : position.yanchor === 'middle' ? height / 2 : 0);
    let score = 0;
    for (let i = 0; i < tt.length; i++) {
      const x = (tt[i] - tt[0]) / (span || 1);
      if (x < left || x > left + width) continue;
      for (const values of [observed, calculated]) {
        const y = mainBottom + (Number(values[i]) || 0) / (maxY * 1.05) * (1 - mainBottom);
        if (y >= bottom && y <= bottom + height) score++;
      }
    }
    if (score < bestScore) { best = position; bestScore = score; }
  }
  return best;
}

async function updateXrdLegendSelection(select) {
  const location = select.value;
  const controls = ['xrd-result-legend-location']
    .map(id => document.getElementById(id)).filter(Boolean);
  controls.forEach(control => { control.value = location; });
  const data = xrdLastResult;
  if (!data) return;
  data.legend_location = location;
  renderXrdPublicationPlot(data);
  const status = document.getElementById('xrd-legend-status');
  if (!data.plot_token) {
    status.textContent = 'Preview updated. Run the fit again to update PNG exports.';
    return;
  }
  status.textContent = 'Updating light and dark PNG exports…';
  controls.forEach(control => { control.disabled = true; });
  try {
    const response = await fetch('/api/xrd/regenerate_plot', {
      method:'POST', headers:{'Content-Type':'application/json'},
      body:JSON.stringify({plot_token:data.plot_token, legend_location:location}),
    });
    const updated = await response.json();
    if (!response.ok || updated.error) throw new Error(updated.error || 'PNG export update failed.');
    if (xrdLastResult !== data) return;
    Object.assign(data, updated, {_plot_revision:Date.now()});
    renderXrdPublicationPlot(data);
    status.textContent = 'Legend updated in the preview and both PNG exports. Fit unchanged.';
  } catch (error) {
    if (xrdLastResult === data) status.textContent = `PNG exports not updated: ${error.message}`;
  } finally {
    controls.forEach(control => { control.disabled = false; });
  }
}
