"""
modules/xrd/cif_cache.py
Persistent disk cache for CIF files.

Stores CIF text files in ~/.catalysis_toolkit_cache/ (or a custom directory
from config.yaml). Cache is keyed by source:id strings, e.g.:
  "cod:1010048"   → COD entry (raw COD CIF)
  "mp:mp-91"      → Materials Project entry (raw MP CIF, OR the
                    fixture-overridden CIF when one is registered in
                    mp_api._LOCAL_FIXTURES)
  "manual:sample" → user-uploaded CIF

SEPARATION RULE (raw vs. prepared):
  Currently this cache holds whichever CIF cached_fetch_mp returned —
  raw MP CIF when no fixture exists, or the fixture content when one
  does.  That is convenient but conflates two roles.  The architecture
  spec (CLAUDE OUTPUTS / Catalysis-Toolkit_Architecture_v1.md, §Step 4)
  recommends separate cache keys for raw / prepared / preview artifacts:
     mp:<id>:raw            — exact MP API response
     mp:<id>:gsas:<vN>      — CIF prepared for GSAS-II refinement
     mp:<id>:preview:<vN>   — preview reflection list
  Future work; the current single-key cache is a known compromise.

Cache never expires — CIF files don't change.
Size is capped at max_size_mb (default 500 MB); oldest files pruned when exceeded.
"""

import os, json, hashlib, time
from pathlib import Path


class CIFCache:
    def __init__(self, cache_dir=None, max_size_mb=500):
        if cache_dir is None:
            cache_dir = os.path.join(Path.home(), '.catalysis_toolkit_cache')
        self.cache_dir   = os.path.expanduser(str(cache_dir))
        self.max_size_mb = max_size_mb
        self.index_path  = os.path.join(self.cache_dir, '_index.json')
        os.makedirs(self.cache_dir, exist_ok=True)
        self._index = self._load_index()

    # ── Index ────────────────────────────────────────────────────────────────

    def _load_index(self):
        if os.path.exists(self.index_path):
            try:
                with open(self.index_path) as f:
                    return json.load(f)
            except Exception:
                pass
        return {}

    def _save_index(self):
        try:
            with open(self.index_path, 'w') as f:
                json.dump(self._index, f)
        except Exception:
            pass

    # ── File naming ──────────────────────────────────────────────────────────

    def _key_to_filename(self, key):
        safe = re.sub(r'[^\w\-]', '_', key)
        return os.path.join(self.cache_dir, f'{safe}.cif')

    # ── Public API ───────────────────────────────────────────────────────────

    def get(self, key):
        """Return cached CIF text for key, or None if not cached."""
        fname = self._key_to_filename(key)
        if key in self._index and os.path.exists(fname):
            try:
                with open(fname, 'r', encoding='utf-8') as f:
                    return f.read()
            except Exception:
                pass
        return None

    def put(self, key, cif_text):
        """Store CIF text under key."""
        fname = self._key_to_filename(key)
        try:
            with open(fname, 'w', encoding='utf-8') as f:
                f.write(cif_text)
            self._index[key] = {
                'file':    fname,
                'size':    len(cif_text.encode('utf-8')),
                'added':   time.time(),
            }
            self._save_index()
            self._maybe_prune()
        except Exception:
            pass

    def has(self, key):
        """Return True if key is in cache and file exists."""
        fname = self._key_to_filename(key)
        return key in self._index and os.path.exists(fname)

    def stats(self):
        """Return dict with entry count and total size in MB."""
        total = sum(v.get('size', 0) for v in self._index.values())
        return {
            'entries':  len(self._index),
            'size_mb':  round(total / (1024**2), 2),
            'max_mb':   self.max_size_mb,
            'dir':      self.cache_dir,
        }

    def clear(self):
        """Delete all cached CIFs and the index."""
        for key, info in list(self._index.items()):
            try:
                os.remove(info['file'])
            except Exception:
                pass
        self._index = {}
        self._save_index()

    # ── Pruning ──────────────────────────────────────────────────────────────

    def _maybe_prune(self):
        total_bytes = sum(v.get('size', 0) for v in self._index.values())
        limit_bytes = self.max_size_mb * 1024 * 1024
        if total_bytes <= limit_bytes:
            return
        # Sort by oldest first
        sorted_keys = sorted(self._index.keys(),
                              key=lambda k: self._index[k].get('added', 0))
        for key in sorted_keys:
            if total_bytes <= limit_bytes * 0.8:  # prune to 80% of limit
                break
            info = self._index.pop(key, {})
            total_bytes -= info.get('size', 0)
            try:
                os.remove(info['file'])
            except Exception:
                pass
        self._save_index()


import re  # needed by _key_to_filename


# ── Module-level singleton ────────────────────────────────────────────────────

_cache_instance = None

def get_cif(key):
    """Return cached CIF text for a key (e.g. 'cod:1010048', 'mp:mp-91'), or None."""
    return get_cache().get(key)


def get_cache(cache_dir=None, max_size_mb=500):
    """Get the shared cache instance (creates on first call)."""
    global _cache_instance
    if _cache_instance is None:
        _cache_instance = CIFCache(cache_dir=cache_dir, max_size_mb=max_size_mb)
    return _cache_instance


def cached_fetch_cod(cod_id, fetch_fn):
    """
    Fetch a COD CIF, using disk cache.
    fetch_fn: callable(cod_id) → parsed_struct_dict (from cod_api.fetch_cif)
    """
    cache = get_cache()
    key   = f'cod:{cod_id}'
    text  = cache.get(key)
    if text:
        from .crystallography import parse_cif
        parsed = parse_cif(text)
        parsed['cod_id']   = str(cod_id)
        parsed['cif_text'] = text
        parsed['cached']   = True
        return parsed
    # Not cached — fetch and store
    result = fetch_fn(cod_id)
    if result.get('cif_text'):
        cache.put(key, result['cif_text'])
    result['cached'] = False
    return result


def cached_fetch_mp(mp_id, api_key, fetch_fn):
    """
    Fetch a Materials Project CIF, using disk cache.
    fetch_fn: callable(mp_id, api_key) → parsed_struct_dict

    Local fixtures (mp_api._LOCAL_FIXTURES) take priority over the disk
    cache.  This ensures that for known-problematic MP entries (e.g.
    mp-2034 W2C, where the round-tripped MP CIF can land as P1/full-cell
    instead of Pbcn), the audited canonical CIF in fixtures/ is always
    what reaches the refinement — regardless of whatever may have been
    cached previously.  The cache is also refreshed with the fixture
    so downstream consumers stay consistent.
    """
    cache = get_cache()
    key   = f'mp:{mp_id}'

    # ── Fixture takes priority over cache ─────────────────────────────
    try:
        from .mp_api import _fixture_cif_for
        fixture_text = _fixture_cif_for(mp_id)
    except Exception:
        fixture_text = None
    if fixture_text:
        from .crystallography import parse_cif
        parsed = parse_cif(fixture_text)
        parsed['mp_id']    = mp_id
        parsed['cod_id']   = mp_id
        parsed['cif_text'] = fixture_text
        parsed['source']   = 'mp'
        parsed['cached']   = False
        # Refresh cache with the fixture so subsequent fetches stay
        # consistent even if this code path ever skips the fixture lookup.
        cache.put(key, fixture_text)
        return parsed

    # No fixture → normal cache flow
    text = cache.get(key)
    if text:
        from .crystallography import parse_cif
        parsed = parse_cif(text)
        parsed['mp_id']    = mp_id
        parsed['cod_id']   = mp_id
        parsed['cif_text'] = text
        parsed['source']   = 'mp'
        parsed['cached']   = True
        return parsed
    result = fetch_fn(mp_id, api_key)
    if result.get('cif_text'):
        cache.put(key, result['cif_text'])
    result['cached'] = False
    return result
