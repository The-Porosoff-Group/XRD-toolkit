#!/usr/bin/env python
"""
Fetch Materials Project CIFs into a local folder for terminal/batch workflows.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def _safe_name(value: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", str(value or "")).strip("_") or "phase"


def _load_api_key(config_path: Path) -> str:
    env_key = os.environ.get("MP_API_KEY") or os.environ.get("MATERIALS_PROJECT_API_KEY")
    if env_key:
        return env_key.strip()
    if not config_path.exists():
        return ""
    try:
        import yaml

        data = yaml.safe_load(_read_text(config_path)) or {}
        for key in ("mp_api_key", "materials_project_api_key", "api_key"):
            if data.get(key):
                return str(data[key]).strip()
        mp = data.get("materials_project") or {}
        if isinstance(mp, dict) and mp.get("api_key"):
            return str(mp["api_key"]).strip()
    except Exception:
        pass
    text = _read_text(config_path)
    m = re.search(r"(?im)^\s*(?:mp_api_key|materials_project_api_key|api_key)\s*:\s*['\"]?([^'\"\s#]+)", text)
    return m.group(1).strip() if m else ""


def fetch_one(mp_id: str, api_key: str, out_dir: Path, overwrite: bool = False) -> Dict[str, Any]:
    from modules.xrd.cif_cache import cached_fetch_mp
    from modules.xrd.mp_api import fetch_cif

    phase = cached_fetch_mp(mp_id, api_key, fetch_cif)
    formula = phase.get("formula") or mp_id
    name = f"{_safe_name(formula)}_{_safe_name(mp_id)}.cif"
    path = out_dir / name
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Use --overwrite to replace it.")
    path.write_text(phase["cif_text"], encoding="utf-8")
    return {
        "mp_id": mp_id,
        "formula": formula,
        "spacegroup": phase.get("spacegroup") or phase.get("spacegroup_name"),
        "spacegroup_number": phase.get("spacegroup_number"),
        "path": str(path),
        "cached": phase.get("cached", False),
    }


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Fetch Materials Project CIFs into a local folder.")
    p.add_argument("--mp-ids", nargs="+", required=True, help="Materials Project ids, e.g. mp-2034")
    p.add_argument("--out-dir", default="cifs", help="Output directory for CIF files")
    p.add_argument("--config", default=str(ROOT / "config.yaml"), help="Config file containing Materials Project API key")
    p.add_argument("--mp-api-key", default="", help="Materials Project API key override")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing CIF files")
    p.add_argument("--manifest", default="manifest.json", help="Manifest filename written inside --out-dir")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    api_key = args.mp_api_key or _load_api_key(Path(args.config))
    rows = []
    for mp_id in args.mp_ids:
        row = fetch_one(mp_id, api_key, out_dir, overwrite=args.overwrite)
        rows.append(row)
        print(f"{row['mp_id']} -> {row['path']}", flush=True)
    manifest_path = out_dir / args.manifest
    manifest_path.write_text(json.dumps({"cifs": rows}, indent=2), encoding="utf-8")
    print(f"Manifest: {manifest_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
