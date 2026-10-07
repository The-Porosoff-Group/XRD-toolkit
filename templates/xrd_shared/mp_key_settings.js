let mpKeyBusy = false;

function updateSourceBadge() {
  const src = document.querySelector('input[name="xrd-source-db"]:checked')?.value || 'both';
  const badge = document.getElementById('mp-key-badge');
  badge.className = window.mpKeySet ? 'pill green' : 'pill yellow';
  badge.textContent = (src === 'mp' || src === 'both')
    ? (window.mpKeySet ? 'MP key saved' : 'MP: no key') : '';
  document.getElementById('mp-key-state').textContent = window.mpKeySet
    ? 'A key is saved. Paste a new key below to replace it.'
    : 'No key saved. Add one to search Materials Project.';
  document.getElementById('mp-api-key').placeholder = window.mpKeySet
    ? 'Paste a replacement API key' : 'Paste your Materials Project API key';
  document.getElementById('mp-key-remove').disabled = mpKeyBusy || !window.mpKeySet;
}

async function manageMpKey(action) {
  if (mpKeyBusy) return;
  const input = document.getElementById('mp-api-key');
  const feedback = document.getElementById('mp-key-feedback');
  const key = input.value.trim();
  if (action === 'save' && !key) {
    feedback.textContent = 'Paste a key to save, or use Remove key to delete the saved key.';
    return;
  }
  mpKeyBusy = true;
  const controls = ['mp-api-key', 'mp-key-save', 'mp-key-test', 'mp-key-remove'];
  controls.forEach(id => { document.getElementById(id).disabled = true; });
  feedback.textContent = action === 'test' ? 'Testing key…' : 'Updating key…';
  try {
    const testing = action === 'test';
    const options = {method: testing ? 'POST' : action === 'remove' ? 'DELETE' : 'PUT'};
    if (action !== 'remove') {
      options.headers = {'Content-Type': 'application/json'};
      options.body = JSON.stringify(key ? {api_key: key} : {});
    }
    const response = await fetch(testing ? '/api/xrd/validate_mp_key' : '/api/xrd/mp_key', options);
    if (response.status === 404 || response.status === 405) {
      throw new Error('Restart the toolkit to enable API key settings, then refresh this page.');
    }
    const data = await response.json();
    if (!response.ok || data.error) throw new Error(data.error || 'Could not update the API key.');
    feedback.textContent = data.message;
    if (!testing) {
      window.mpKeySet = data.mp_key_set;
      input.value = '';
      updateSourceBadge();
    }
  } catch (error) {
    feedback.textContent = error.message;
  } finally {
    mpKeyBusy = false;
    controls.forEach(id => { document.getElementById(id).disabled = false; });
    updateSourceBadge();
  }
}
