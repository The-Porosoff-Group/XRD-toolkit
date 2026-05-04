#!/usr/bin/env python
"""
Batch XRD refinement runner.

Runs the same modules.xrd.run(...) backend used by the GUI, but reads phases
and refinement controls from a GUI preset, a CLI recipe JSON file, CIF files,
and/or Materials Project ids.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    return data


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def _safe_name(value: str) -> str:
    safe = re.sub(r"[^\w.\-]+", "_", str(value or "").strip())
    return safe.strip("_") or "sample"


def _norm(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes", "on"}


def _as_float(value: Any, default: Optional[float] = None) -> Optional[float]:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


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


def _builtin_presets() -> List[Dict[str, Any]]:
    return [{
        "id": "builtin-wc-w2c-synergy-s",
        "name": "WC/W2C Synergy-S production",
        "recipe_key": "wc_w2c_synergy_s",
        "recipe_only": True,
        "controls": {
            "wavelength": "1.54056",
            "tt_min": "20",
            "tt_max": "60",
            "n_bg_coeffs": "auto",
            "instrument": "synergy_s",
            "checkboxes": {
                "xrd-calibration-mode": False,
                "xrd-verification-mode": True,
                "xrd-verify-cell": True,
                "xrd-phase-isolation": True,
                "xrd-zero-not-disp": True,
                "xrd-refine-x": True,
                "xrd-fix-y": False,
                "xrd-y-nonneg": True,
                "xrd-refine-uiso": False,
                "xrd-refine-xyz": False,
            },
        },
        "phase_option_rules": [
            {
                "match": "WC",
                "refine_size": False,
                "refine_mustrain": False,
                "po_mode": "fixed",
                "po_axis": [0, 0, 1],
                "po_value": 0.905,
                "uniform_cell": False,
            },
            {
                "match": "W2C",
                "refine_size": False,
                "refine_mustrain": False,
                "po_mode": "off",
                "po_axis": [0, 0, 1],
                "po_value": 1.0,
                "uniform_cell": False,
            },
        ],
    }]


def _load_presets(path: Path) -> List[Dict[str, Any]]:
    presets = list(_builtin_presets())
    if not path.exists():
        return presets
    data = _load_json(path)
    raw = data.get("presets", data)
    if isinstance(raw, list):
        presets.extend(p for p in raw if isinstance(p, dict))
    return presets


def _select_preset(path: Path, selector: str) -> Dict[str, Any]:
    wanted = _norm(selector)
    for preset in _load_presets(path):
        if wanted in {_norm(preset.get("id")), _norm(preset.get("name"))}:
            return dict(preset)
    names = ", ".join(p.get("name", p.get("id", "?")) for p in _load_presets(path))
    raise ValueError(f"Preset '{selector}' not found in {path}. Available: {names}")


def _phase_from_cif(path: Path) -> Dict[str, Any]:
    from modules.xrd.crystallography import parse_cif

    text = _read_text(path)
    parsed = parse_cif(text)
    name = path.stem
    parsed.update({
        "name": parsed.get("formula") or name,
        "formula": parsed.get("formula") or name,
        "cod_id": "manual",
        "source": "manual",
        "cif_text": text,
        "cif_path": str(path),
    })
    return parsed


def _fetch_mp_phase(mp_id: str, api_key: str) -> Dict[str, Any]:
    from modules.xrd.cif_cache import cached_fetch_mp
    from modules.xrd.mp_api import fetch_cif

    return cached_fetch_mp(mp_id, api_key, fetch_cif)


def _phase_keys(phase: Dict[str, Any]) -> set:
    return {
        _norm(phase.get("name")),
        _norm(phase.get("formula")),
        _norm(phase.get("mp_id")),
        _norm(phase.get("cod_id")),
        _norm(Path(str(phase.get("cif_path", ""))).stem),
    }


def _load_cif_dir(cif_dir: Optional[Path]) -> List[Dict[str, Any]]:
    if not cif_dir:
        return []
    if not cif_dir.exists():
        raise FileNotFoundError(f"CIF directory does not exist: {cif_dir}")
    return [_phase_from_cif(p) for p in sorted(cif_dir.glob("*.cif"))]


def _resolve_phases(
    recipe: Dict[str, Any],
    recipe_path: Optional[Path],
    cif_dir: Optional[Path],
    cif_files: Iterable[Path],
    mp_ids: Iterable[str],
    api_key: str,
) -> List[Dict[str, Any]]:
    base_dir = recipe_path.parent if recipe_path else ROOT
    phases: List[Dict[str, Any]] = []

    for raw in recipe.get("phases") or []:
        if not isinstance(raw, dict):
            continue
        ph = dict(raw)
        cif_path = ph.get("cif_path") or ph.get("cif_file")
        if cif_path and not ph.get("cif_text"):
            path = Path(cif_path)
            if not path.is_absolute():
                path = base_dir / path
            parsed = _phase_from_cif(path)
            parsed.update({k: v for k, v in ph.items() if v not in (None, "")})
            ph = parsed
        if (ph.get("mp_id") or str(ph.get("cod_id", "")).startswith("mp-")) and not ph.get("cif_text"):
            mp_id = str(ph.get("mp_id") or ph.get("cod_id"))
            fetched = _fetch_mp_phase(mp_id, api_key)
            fetched.update({k: v for k, v in ph.items() if v not in (None, "")})
            ph = fetched
        phases.append(ph)

    loaded_cifs = _load_cif_dir(cif_dir)
    loaded_cifs.extend(_phase_from_cif(Path(p)) for p in cif_files)

    if not phases:
        phases.extend(loaded_cifs)
    else:
        unused = []
        for cif_phase in loaded_cifs:
            matched = False
            cif_keys = _phase_keys(cif_phase)
            for i, ph in enumerate(phases):
                if _phase_keys(ph) & cif_keys:
                    merged = dict(cif_phase)
                    merged.update({k: v for k, v in ph.items() if v not in (None, "")})
                    if not merged.get("cif_text"):
                        merged["cif_text"] = cif_phase.get("cif_text")
                    phases[i] = merged
                    matched = True
                    break
            if not matched:
                unused.append(cif_phase)
        phases.extend(unused)

    for mp_id in mp_ids:
        phases.append(_fetch_mp_phase(mp_id, api_key))

    if not phases:
        raise ValueError("No phases found. Provide --preset with phases, --recipe, --cif-dir, --cif, or --mp-ids.")
    return phases


def _axis(value: Any) -> List[int]:
    if isinstance(value, (list, tuple)):
        vals = value
    else:
        vals = str(value or "0 0 1").replace(",", " ").split()
    out = []
    for v in vals[:3]:
        try:
            out.append(int(float(v)))
        except (TypeError, ValueError):
            out.append(0)
    while len(out) < 3:
        out.append(0)
    return out


def _default_phase_option() -> Dict[str, Any]:
    return {
        "refine_size": False,
        "refine_mustrain": False,
        "po_mode": "off",
        "po_axis": [0, 0, 1],
        "po_value": 1.0,
        "uniform_cell": False,
    }


def _normalize_phase_options(recipe: Dict[str, Any], phases: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    raw_options = recipe.get("phase_options")
    if isinstance(raw_options, list) and raw_options:
        options = []
        for opt in raw_options[:len(phases)]:
            merged = _default_phase_option()
            if isinstance(opt, dict):
                merged.update(opt)
            merged["po_axis"] = _axis(merged.get("po_axis"))
            merged["po_value"] = _as_float(merged.get("po_value"), 1.0)
            options.append(merged)
        while len(options) < len(phases):
            options.append(_default_phase_option())
        return options

    options = [_default_phase_option() for _ in phases]
    for rule in recipe.get("phase_option_rules") or []:
        if not isinstance(rule, dict):
            continue
        match = _norm(rule.get("match") or rule.get("formula") or rule.get("name") or rule.get("mp_id"))
        if not match:
            continue
        for i, phase in enumerate(phases):
            if match in _phase_keys(phase):
                merged = dict(options[i])
                merged.update({k: v for k, v in rule.items() if k not in {"match", "formula", "name"}})
                merged["po_axis"] = _axis(merged.get("po_axis"))
                merged["po_value"] = _as_float(merged.get("po_value"), 1.0)
                options[i] = merged
    return options


def _controls_to_params(recipe: Dict[str, Any], phases: List[Dict[str, Any]], instprm_file: Optional[Path]) -> Dict[str, Any]:
    controls = recipe.get("controls") or {}
    cb = controls.get("checkboxes") or {}

    params = {
        "phases": phases,
        "wavelength": _as_float(controls.get("wavelength"), 1.54056),
        "wavelength_label": controls.get("wavelength_label") or controls.get("wavelength_source") or "",
        "tt_min": _as_float(controls.get("tt_min"), None),
        "tt_max": _as_float(controls.get("tt_max"), None),
        "n_bg_coeffs": controls.get("n_bg_coeffs", "auto"),
        "max_outer": int(_as_float(controls.get("max_outer"), 10) or 10),
        "method": recipe.get("method") or controls.get("method") or "gsas2",
        "instprm_file": str(instprm_file) if instprm_file else controls.get("instprm_file"),
        "instrument": controls.get("instrument", "auto"),
        "verification_mode": _as_bool(cb.get("xrd-verification-mode")),
        "verify_refine_cell": _as_bool(cb.get("xrd-verify-cell")),
        "phase_isolation": _as_bool(cb.get("xrd-phase-isolation")),
        "verify_refine_po": _as_bool(cb.get("xrd-verify-refine-po")),
        "verify_use_zero_not_displace": _as_bool(cb.get("xrd-zero-not-disp")),
        "verify_cell_uniform_w2c": _as_bool(cb.get("xrd-uniform-w2c-cell")),
        "verify_refine_x": _as_bool(cb.get("xrd-refine-x")),
        "verify_fix_y": _as_bool(cb.get("xrd-fix-y")),
        "verify_y_fixed_value": _as_float(controls.get("fix_y_value"), None),
        "verify_y_nonnegative": _as_bool(cb.get("xrd-y-nonneg")),
        "verify_refine_uiso": _as_bool(cb.get("xrd-refine-uiso")),
        "verify_refine_size": _as_bool(cb.get("xrd-refine-size")),
        "refine_xyz": _as_bool(cb.get("xrd-refine-xyz")),
        "use_gsas_ref_ticks": _as_bool(cb.get("xrd-use-gsas-ref-ticks")),
        "verify_fix_po": _as_bool(cb.get("xrd-fix-po")),
        "verify_po_fixed_value": _as_float(controls.get("po_fixed_value"), None),
        "phase_options": _normalize_phase_options(recipe, phases),
    }
    if params["tt_min"] is None:
        params.pop("tt_min")
    if params["tt_max"] is None:
        params.pop("tt_max")
    return params


def _expand_patterns(patterns: Iterable[str]) -> List[Path]:
    files: List[Path] = []
    for pattern in patterns:
        hits = glob.glob(pattern)
        if hits:
            files.extend(Path(h) for h in hits)
        else:
            files.append(Path(pattern))
    files = [p for p in files if p.exists() and p.is_file()]
    if not files:
        raise FileNotFoundError("No XRD pattern files matched.")
    return sorted(files)


def _compact_phase(sample: str, stats: Dict[str, Any], phase: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "sample": sample,
        "Rwp": stats.get("Rwp"),
        "GoF": stats.get("GoF"),
        "phase": phase.get("name") or phase.get("formula"),
        "formula": phase.get("formula"),
        "spacegroup": phase.get("spacegroup"),
        "weight_fraction_pct": phase.get("weight_fraction_%"),
        "weight_fraction_err_pct": phase.get("weight_fraction_err_%"),
        "weight_fraction_err_source": phase.get("weight_fraction_err_source"),
        "preferred_orientation_value": phase.get("preferred_orientation_value"),
        "preferred_orientation_axis": " ".join(map(str, phase.get("preferred_orientation_axis") or [])),
        "preferred_orientation_mode": phase.get("preferred_orientation_mode"),
        "crystallite_size_nm": phase.get("crystallite_size_nm"),
        "crystallite_size_source": phase.get("crystallite_size_source"),
        "fwhm_deg": phase.get("fwhm_deg"),
        "fwhm_reference_two_theta": phase.get("fwhm_reference_two_theta"),
        "delta_a_pct": phase.get("delta_a_pct"),
        "delta_b_pct": phase.get("delta_b_pct"),
        "delta_c_pct": phase.get("delta_c_pct"),
        "delta_volume_pct": phase.get("delta_volume_pct"),
    }


def _write_outputs(output_dir: Path, rows: List[Dict[str, Any]], results: List[Dict[str, Any]]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "batch_phase_summary.csv"
    if rows:
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    json_path = output_dir / "batch_summary.json"
    with json_path.open("w", encoding="utf-8") as f:
        json.dump({"created_at": datetime.now().isoformat(timespec="seconds"), "results": results}, f, indent=2)


def run_batch(args: argparse.Namespace) -> int:
    from modules import xrd as xrd_module

    recipe_path = Path(args.recipe).resolve() if args.recipe else None
    if recipe_path:
        recipe = _load_json(recipe_path)
    elif args.preset:
        recipe = _select_preset(Path(args.presets_file), args.preset)
    else:
        recipe = _select_preset(Path(args.presets_file), "builtin-wc-w2c-synergy-s")

    api_key = args.mp_api_key or _load_api_key(Path(args.config))
    phases = _resolve_phases(
        recipe=recipe,
        recipe_path=recipe_path,
        cif_dir=Path(args.cif_dir).resolve() if args.cif_dir else None,
        cif_files=[Path(p).resolve() for p in args.cif],
        mp_ids=args.mp_ids,
        api_key=api_key,
    )
    params = _controls_to_params(recipe, phases, Path(args.instprm).resolve() if args.instprm else None)
    patterns = _expand_patterns(args.patterns)

    output_root = Path(args.out).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    resolved_recipe = {
        "recipe": recipe,
        "resolved_phase_names": [p.get("name") or p.get("formula") for p in phases],
        "params": {k: v for k, v in params.items() if k != "phases"},
    }
    with (output_root / "resolved_batch_recipe.json").open("w", encoding="utf-8") as f:
        json.dump(resolved_recipe, f, indent=2)

    all_rows: List[Dict[str, Any]] = []
    batch_results: List[Dict[str, Any]] = []
    for pattern in patterns:
        sample = args.sample_id or pattern.stem
        sample_dir = output_root / _safe_name(sample)
        if len(patterns) > 1:
            sample_dir = output_root / _safe_name(pattern.stem)
            sample = pattern.stem
        print(f"\n=== XRD batch: {pattern.name} ===", flush=True)
        result = xrd_module.run(
            filepath=str(pattern),
            output_dir=str(sample_dir),
            metadata={"sample_id": sample, "notes": f"Batch recipe: {recipe.get('name', recipe.get('id', 'recipe'))}"},
            params=params,
        )
        stats = result.get("statistics", {})
        rows = [_compact_phase(sample, stats, ph) for ph in result.get("phase_results", [])]
        all_rows.extend(rows)
        compact = {
            "sample": sample,
            "pattern": str(pattern),
            "output_dir": str(sample_dir),
            "statistics": stats,
            "phase_results": rows,
            "plot_path": result.get("plot_path"),
            "summary_path": result.get("summary_path"),
        }
        batch_results.append(compact)
        with (sample_dir / "summary.json").open("w", encoding="utf-8") as f:
            json.dump(compact, f, indent=2)
        print(f"Rwp={stats.get('Rwp')}  GoF={stats.get('GoF')}  output={sample_dir}", flush=True)

    _write_outputs(output_root, all_rows, batch_results)
    print(f"\nBatch complete: {output_root}", flush=True)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run XRD refinements from the terminal using GUI presets, recipes, CIFs, or MP ids.")
    p.add_argument("--patterns", nargs="+", required=True, help="XRD pattern file(s) or glob(s), e.g. data/*.xy")
    p.add_argument("--out", default="results/xrd_batch", help="Batch output directory")
    p.add_argument("--preset", help="Preset name or id from xrd_refinement_presets.json; built-in WC/W2C is available")
    p.add_argument("--presets-file", default=str(ROOT / "xrd_refinement_presets.json"), help="GUI preset JSON path")
    p.add_argument("--recipe", help="Recipe JSON path. Overrides --preset.")
    p.add_argument("--cif-dir", help="Directory of .cif files to use as phases")
    p.add_argument("--cif", action="append", default=[], help="One CIF file to add. Can be repeated.")
    p.add_argument("--mp-ids", nargs="*", default=[], help="Materials Project ids to fetch as phases, e.g. mp-2034")
    p.add_argument("--config", default=str(ROOT / "config.yaml"), help="Config file containing Materials Project API key")
    p.add_argument("--mp-api-key", default="", help="Materials Project API key override")
    p.add_argument("--instprm", help="Optional GSAS-II .instprm file")
    p.add_argument("--sample-id", help="Sample id for a single pattern run")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return run_batch(args)


if __name__ == "__main__":
    raise SystemExit(main())
