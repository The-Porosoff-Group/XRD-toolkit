"""Audit XRD CIF fixtures against fixture_manifest.json.

This is intentionally conservative. It does not regenerate CIFs and it does
not "fix" headers. Its job is to make it obvious which fixtures are canonical
normal-import CIFs and which are raw/regression source files that need special
handling.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "fixtures"
MANIFEST = FIXTURES / "fixture_manifest.json"


def _load_parser():
    sys.path.insert(0, str(ROOT))
    from modules.xrd.crystallography import parse_cif

    return parse_cif


def _norm_sg(value):
    return str(value or "").replace(" ", "").replace("_", "").lower()


def audit(strict: bool = False) -> int:
    parse_cif = _load_parser()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    entries = manifest.get("fixtures", {})
    failures = []
    warnings = []

    for path in sorted(FIXTURES.glob("*.cif")):
        rec = entries.get(path.name)
        if not rec:
            warnings.append(f"{path.name}: missing from fixture_manifest.json")
            continue

        parsed = parse_cif(path.read_text(encoding="utf-8", errors="replace"))
        parsed_sg = int(parsed.get("spacegroup_number") or 1)
        expected_sg = int(rec.get("spacegroup_number") or parsed_sg)
        intended = set(rec.get("intended_use") or [])
        normal_safe = rec.get("normal_import_safe", True)

        if "normal_import" in intended and not normal_safe:
            failures.append(f"{path.name}: marked normal_import but normal_import_safe=false")

        if "normal_import" in intended and parsed_sg != expected_sg:
            failures.append(
                f"{path.name}: normal_import fixture parses as SG {parsed_sg}, "
                f"expected SG {expected_sg}"
            )
        elif parsed_sg != expected_sg:
            warnings.append(
                f"{path.name}: parses as SG {parsed_sg}, manifest intent is SG "
                f"{expected_sg} ({rec.get('cell_setting', 'unknown')})"
            )

        parsed_name = parsed.get("spacegroup") or parsed.get("spacegroup_name")
        if "normal_import" in intended and parsed_name and rec.get("spacegroup"):
            if _norm_sg(parsed_name) != _norm_sg(rec.get("spacegroup")):
                failures.append(
                    f"{path.name}: H-M symbol '{parsed_name}' does not match "
                    f"manifest '{rec.get('spacegroup')}'"
                )

    for name in sorted(entries):
        if not (FIXTURES / name).exists():
            warnings.append(f"{name}: listed in manifest but file is missing")

    if warnings:
        print("Fixture audit warnings:")
        for msg in warnings:
            print(f"  - {msg}")
    if failures:
        print("Fixture audit failures:")
        for msg in failures:
            print(f"  - {msg}")
        return 1
    if strict and warnings:
        return 1
    print("Fixture audit completed.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict", action="store_true", help="Treat warnings as failures")
    args = parser.parse_args()
    return audit(strict=args.strict)


if __name__ == "__main__":
    raise SystemExit(main())
