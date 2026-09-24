const XRD_INSTRUMENT_IDS = ['smartlab', 'synergy_s', 'benchtop_cu', 'none'];

function isXrdCalibration() {
  return document.getElementById('xrd-instrument').value === 'none';
}

function xrdInstrumentChanged(clearUpload = true) {
  const calibrating = isXrdCalibration();
  const upload = document.getElementById('xrd-instprm-file');
  if (clearUpload) upload.value = '';
  upload.disabled = calibrating;
  document.getElementById('xrd-calibration-help').style.display = calibrating ? '' : 'none';
  document.getElementById('xrd-instrument-status').textContent = '';
  document.getElementById('btn-xrd-gsas2').textContent = calibrating ? 'Run calibration' : 'GSAS-II Refinement';
  document.getElementById('xrd-instprm-override').style.display = calibrating ? 'none' : '';
  document.getElementById('xrd-calibration-phase-hint').style.display = calibrating ? '' : 'none';
}

function restoreXrdInstrument(controls = {}) {
  const selected = controls.instrument;
  const legacyCalibration = controls.checkboxes?.['xrd-calibration-mode'];
  const supported = XRD_INSTRUMENT_IDS.includes(selected);
  document.getElementById('xrd-instrument').value = legacyCalibration ? 'none' : supported ? selected : 'none';
  document.getElementById('xrd-calibration-geometry').value =
    selected === 'generic_capillary' || (legacyCalibration && selected === 'synergy_s') ? 'capillary' :
    legacyCalibration && ['smartlab', 'benchtop_cu'].includes(selected) ? 'bragg_brentano' :
    controls.instrument_geometry || 'bragg_brentano';
  xrdInstrumentChanged();
  if (selected && !supported && !legacyCalibration) {
    document.getElementById('xrd-instrument-status').textContent =
      'This saved preset used an older profile. Choose one of the three instruments for sample fitting, or None / calibration for a standard.';
  }
}

function appendXrdInstrumentSettings(fd) {
  const selected = document.getElementById('xrd-instrument').value;
  if (!XRD_INSTRUMENT_IDS.includes(selected)) throw new Error('Select an instrument or None / calibration.');
  fd.append('instrument', selected);
  if (isXrdCalibration()) {
    fd.append('calibration_mode', 'true');
    fd.append('instrument_geometry', document.getElementById('xrd-calibration-geometry').value);
    fd.append('calibration_standard', 'Si640g');
  } else {
    const upload = document.getElementById('xrd-instprm-file').files[0];
    if (upload) fd.append('instprm_file', upload);
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
    : 'Candidate created, but parameter checks failed. Review the report and rerun before using this file.';
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
}

document.addEventListener('DOMContentLoaded', () => xrdInstrumentChanged(false));
