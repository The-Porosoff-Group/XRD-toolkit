let xrdInstrumentProfiles = [];

async function refreshXrdInstruments(selected) {
  const select = document.getElementById('xrd-instrument');
  const previous = selected || select.value || 'generic_flat_plate';
  try {
    const response = await fetch('/api/xrd/instruments');
    if (response.status === 404) throw new Error('Restart the toolkit to enable local instrument profiles.');
    const data = await response.json();
    if (!response.ok) throw new Error(data.error || 'Could not load instruments.');
    xrdInstrumentProfiles = data.instruments;
    const options = data.instruments.map(profile => {
      const option = document.createElement('option');
      option.value = profile.id;
      option.textContent = profile.label + (profile.local ? ' — local' : '');
      return option;
    });
    for (const [value, label] of [['upload', 'Upload .instprm…'], ['auto', 'Auto-detect from scan name']]) {
      const option = document.createElement('option');
      option.value = value; option.textContent = label; options.push(option);
    }
    select.replaceChildren(...options);
    select.value = options.some(option => option.value === previous) ? previous : data.default;
    xrdInstrumentChanged(false);
  } catch (error) {
    document.getElementById('xrd-instrument-status').textContent = error.message;
  }
}

function xrdInstrumentChanged(clearUpload = true) {
  const value = document.getElementById('xrd-instrument').value;
  const uploading = value === 'upload';
  document.getElementById('xrd-upload-profile').style.display = uploading ? '' : 'none';
  document.getElementById('xrd-save-upload').style.display = uploading ? '' : 'none';
  if (clearUpload && !uploading) document.getElementById('xrd-instprm-file').value = '';
  const profile = xrdInstrumentProfiles.find(item => item.id === value);
  if (clearUpload && Number.isFinite(profile?.polariz)) {
    document.getElementById('xrd-calibration-polariz').value = profile.polariz;
  }
  let description = profile?.notes || (uploading
    ? 'Choose the geometry and upload your own GSAS-II profile.'
    : 'Unrecognized scans use the generic flat-plate profile, never a named instrument calibration.');
  if (profile?.calibration_range) description += ` Calibrated over ${profile.calibration_range.join('–')}° 2θ.`;
  if (profile?.wavelength) description += ` File wavelength: ${profile.wavelength.toFixed(6)} Å.`;
  document.getElementById('xrd-instrument-description').textContent = description;
}

function xrdCalibrationChanged() {
  const calibrating = document.getElementById('xrd-calibration-mode').checked;
  document.getElementById('xrd-calibration-help').style.display = calibrating ? '' : 'none';
  document.getElementById('xrd-save-upload').disabled = calibrating;
}

function appendXrdInstrumentSettings(fd) {
  const selected = document.getElementById('xrd-instrument').value;
  if (!selected) throw new Error('Select an available instrument or a generic geometry.');
  const calibrating = document.getElementById('xrd-calibration-mode').checked;
  const upload = document.getElementById('xrd-instprm-file').files[0];
  if (selected === 'upload' && !upload && !calibrating) throw new Error('Choose a .instprm file to upload.');
  fd.append('instrument', selected);
  fd.append('instrument_geometry', document.getElementById('xrd-upload-geometry').value);
  fd.append('spectrum', document.getElementById('xrd-spectrum').value);
  fd.append('calibration_standard', 'Si640g');
  fd.append('calibration_name', document.getElementById('xrd-local-instrument-name').value.trim());
  fd.append('calibration_polariz', document.getElementById('xrd-calibration-polariz').value);
  if (selected === 'upload' && upload && !calibrating) fd.append('instprm_file', upload);
}

async function saveXrdLocalInstrument(fromCalibration) {
  const status = document.getElementById(fromCalibration ? 'xrd-calibration-save-status' : 'xrd-instrument-status');
  const button = document.getElementById(fromCalibration ? 'xrd-save-calibration' : 'xrd-save-upload');
  button.disabled = true;
  try {
    const fd = new FormData();
    const input = document.getElementById(fromCalibration ? 'xrd-result-instrument-name' : 'xrd-local-instrument-name');
    const label = input.value.trim();
    if (!label) throw new Error('Enter a name for this local instrument.');
    fd.append('label', label);
    if (fromCalibration) {
      if (!xrdLastResult?.calibration_token) throw new Error('Run calibration first.');
      fd.append('calibration_token', xrdLastResult.calibration_token);
    } else {
      const file = document.getElementById('xrd-instprm-file').files[0];
      if (!file) throw new Error('Choose a .instprm file first.');
      fd.append('instprm_file', file);
      fd.append('geometry', document.getElementById('xrd-upload-geometry').value);
    }
    const response = await fetch('/api/xrd/instruments', {method:'POST', body:fd});
    const data = await response.json();
    if (!response.ok) throw new Error(data.error || 'Could not save instrument.');
    await refreshXrdInstruments(data.instrument);
    document.getElementById('xrd-local-instrument-name').value = label;
    document.getElementById('xrd-instprm-file').value = '';
    if (fromCalibration) {
      document.getElementById('xrd-calibration-mode').checked = false;
      xrdCalibrationChanged();
    }
    status.textContent = `${label} saved locally and selected. New sample fits will load its .instprm automatically.`;
  } catch (error) {
    status.textContent = error.message;
  } finally {
    button.disabled = false;
  }
}

function renderXrdCalibrationResult(data) {
  const panel = document.getElementById('xrd-calibration-result');
  panel.replaceChildren();
  panel.style.display = data.calibration_token ? '' : 'none';
  for (const id of ['xrd-baseline-comparison', 'xrd-reusable-fit', 'xrd-phase-results', 'xrd-pymatgen-pill']) {
    document.getElementById(id).style.display = data.calibration_token ? 'none' : '';
  }
  if (!data.calibration_token) return;
  const validation = data.calibration_validation || {};
  const title = document.createElement('p');
  title.textContent = validation.passed
    ? 'Candidate created. Parameter checks passed; review the fit before using this profile.'
    : 'Candidate created, but parameter checks failed. Review the report and rerun before saving.';
  panel.append(title);
  const details = document.createElement('p');
  details.textContent = [...(validation.reasons || []), ...(validation.warnings || [])].join(' ');
  panel.append(details);
  const standard = document.createElement('p');
  standard.textContent = 'NIST Si 640g · certified a = 5.431109 Å (fixed) · ' +
    (data.geometry === 'bragg_brentano' ? 'flat plate / Bragg-Brentano' : 'capillary / transmission');
  panel.append(standard);
  const parameters = document.createElement('p');
  parameters.textContent = Object.entries(data.phase_results?.[0]?.instprm_params || {})
    .map(([key, value]) => `${key} = ${Number(value).toFixed(5)}`).join(' · ');
  panel.append(parameters);
  for (const [path, label] of [[data.candidate_instprm_path, 'Download candidate .instprm'],
      [data.calibration_report_txt, 'Download calibration report'],
      [data.calibration_report_json, 'Download detailed report (JSON)']]) {
    if (!path) continue;
    const button = document.createElement('button');
    button.type = 'button'; button.className = 'btn-sm'; button.textContent = label;
    button.style.margin = '0 6px 8px 0'; button.onclick = () => downloadFile(path);
    panel.append(button);
  }
  const label = document.createElement('label');
  label.htmlFor = 'xrd-result-instrument-name'; label.textContent = 'Save this candidate under a local instrument name';
  label.style.display = 'block'; panel.append(label);
  const input = document.createElement('input');
  input.id = 'xrd-result-instrument-name'; input.maxLength = 100;
  input.value = data.calibration_name || ''; panel.append(input);
  const save = document.createElement('button');
  save.id = 'xrd-save-calibration'; save.className = 'btn-sm green'; save.type = 'button';
  save.textContent = 'Save as local instrument'; save.disabled = !validation.passed;
  save.onclick = () => saveXrdLocalInstrument(true); panel.append(save);
  const status = document.createElement('p'); status.id = 'xrd-calibration-save-status';
  status.setAttribute('role', 'status'); panel.append(status);
}

document.addEventListener('DOMContentLoaded', () => { refreshXrdInstruments(); xrdCalibrationChanged(); });
