"""Publication-facing labels and reflection metadata for XRD exports.

This module is deliberately display-only.  It formats values that have already
been refined and attaches Miller-index labels to the tick positions used in the
plots and workbook exports; it never changes refinement inputs or calculations.
"""

from __future__ import annotations

from datetime import date, datetime
import os
import re
import unicodedata


_SUBSCRIPT = str.maketrans({
    "0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄",
    "5": "₅", "6": "₆", "7": "₇", "8": "₈", "9": "₉",
    "+": "₊", "-": "₋", "=": "₌", "(": "₍", ")": "₎",
    "x": "ₓ",
})


def format_chemical_formula(formula):
    """Return a plain-Unicode chemical formula with true subscripts."""
    text = str(formula or "").strip().replace("_", "")
    if not text:
        return ""

    # Formula fields contain stoichiometry rather than sample names.  Translating
    # all stoichiometric digits keeps labels portable across Matplotlib, Plotly,
    # PNG metadata, and Excel without relying on a TeX installation.
    text = text.translate(_SUBSCRIPT)
    return text.replace(" ", "")


def format_space_group(space_group):
    """Format Hermann-Mauguin symbols without raw underscore/minus notation."""
    text = str(space_group or "").strip()
    if not text:
        return ""

    # Screw-axis subscripts are commonly supplied as P6_3/mmc.
    text = re.sub(
        r"_([0-9]+)",
        lambda match: match.group(1).translate(_SUBSCRIPT),
        text,
    )
    # A leading minus before an index denotes an overbar, not punctuation.
    text = re.sub(r"-([0-9])", lambda match: match.group(1) + "\u0305", text)
    return re.sub(r"\s+", " ", text).strip()


def clean_descriptive_text(value, fallback="Phase"):
    """Expand filename-style separators for human-facing labels."""
    text = str(value or "").strip()
    text = re.sub(r"[_-]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text or fallback


def phase_display_label(phase, index=None):
    """Return a clear phase label with formula and explicit space-group text."""
    formula = format_chemical_formula(phase.get("formula", ""))
    if formula:
        base = formula
    else:
        fallback = f"Phase {index + 1}" if index is not None else "Phase"
        base = clean_descriptive_text(phase.get("name", ""), fallback=fallback)

    space_group = format_space_group(phase.get("spacegroup", ""))
    if not space_group:
        number = phase.get("spacegroup_number")
        space_group = f"number {number}" if number not in (None, "") else ""
    if space_group:
        return f"{base} (space group {space_group})"
    return base


def phase_legend_label(phase, index=None, *, show_weight=True):
    """Return the publication legend text for one fitted phase."""
    label = phase_tick_label(phase, index=index)
    weight = phase.get("weight_fraction_%")
    uncertainty = phase.get("weight_fraction_err_%")
    if show_weight and weight not in (None, ""):
        value = f"{weight}"
        if uncertainty not in (None, ""):
            value += f" ± {uncertainty}"
        label += f", {value} wt. %"
    return label


def phase_tick_label(phase, index=None):
    """Return a compact formula-and-space-group reflection-row label."""
    formula = format_chemical_formula(phase.get("formula", ""))
    if not formula:
        fallback = f"Fitted phase {index + 1}" if index is not None else "Fitted phase"
        formula = clean_descriptive_text(phase.get("name", ""), fallback)
    space_group = format_space_group(phase.get("spacegroup", ""))
    return f"{formula} ({space_group})" if space_group else formula


def format_wavelength_label(value):
    """Normalize wavelength labels such as Cu Kalpha2 to Cu Kα₂."""
    text = clean_descriptive_text(value, fallback="")
    text = re.sub(r"K\s*(?:alpha|α)\s*([12])?", lambda match: (
        "Kα" + (match.group(1).translate(_SUBSCRIPT) if match.group(1) else "")
    ), text, flags=re.IGNORECASE)
    return text


def format_hkl(hkl):
    """Return conventional spaced Miller indices, including overbars."""
    if isinstance(hkl, str):
        values = re.findall(r"-?\d+", hkl)
    else:
        values = list(hkl or [])
    if not values:
        return ""

    rendered = []
    for value in values[:4]:
        try:
            number = int(value)
        except (TypeError, ValueError):
            rendered.append(str(value))
            continue
        digits = str(abs(number))
        if number < 0:
            digits = "".join(char + "\u0305" for char in digits)
        rendered.append(digits)
    return f"({' '.join(rendered)})"


def _normalize_tick_reflections(raw):
    normalized = []
    for item in raw or []:
        if isinstance(item, dict):
            two_theta = item.get("two_theta", item.get("position"))
            hkl = item.get("hkl")
        elif isinstance(item, (list, tuple)) and len(item) >= 2:
            two_theta, hkl = item[0], item[1]
        else:
            continue
        try:
            two_theta = round(float(two_theta), 3)
        except (TypeError, ValueError):
            continue
        label = format_hkl(hkl)
        if label:
            normalized.append({
                "two_theta": two_theta,
                "hkl": list(hkl) if not isinstance(hkl, str) else hkl,
                "label": label,
            })
    return normalized


def reflection_labels_for_phase(phase, wavelength, tt_min, tt_max):
    """Return the plotted tick positions with their matching Miller indices.

    Backends now provide ``tick_reflections`` directly.  The calculated fallback
    mirrors the historical Excel exporter so older saved results remain usable.
    """
    supplied = _normalize_tick_reflections(phase.get("tick_reflections", []))
    if supplied:
        return supplied

    try:
        from .crystallography import generate_reflections, parse_cif

        system = (phase.get("system") or "triclinic").lower()
        space_group_number = phase.get("spacegroup_number", 1)
        sites = None
        cif_text = phase.get("cif_text", "")
        if cif_text:
            try:
                sites = parse_cif(cif_text).get("sites") or None
            except Exception:
                sites = None

        site_policy = "auto"
        try:
            from .gsasii_backend import _cif_policy
            if _cif_policy(phase) == "mp_w2c_pbcn_compat":
                site_policy = "legacy_direct_sites"
        except Exception:
            pass

        generated = generate_reflections(
            phase.get("a", 1), phase.get("b", 1), phase.get("c", 1),
            phase.get("alpha", 90), phase.get("beta", 90),
            phase.get("gamma", 90), system, space_group_number,
            wavelength, tt_min, tt_max, hkl_max=12,
            sites=sites, site_policy=site_policy,
        )
        filtered_ticks = sorted(float(value) for value in (
            phase.get("tick_positions", []) or []))
        chosen = []
        if filtered_ticks:
            for tick in filtered_ticks:
                nearby = [ref for ref in generated if abs(float(ref[0]) - tick) < 0.02]
                if nearby:
                    ref = min(nearby, key=lambda candidate: abs(float(candidate[0]) - tick))
                    chosen.append((tick, ref[2]))
        else:
            chosen = [(float(ref[0]), ref[2]) for ref in generated]
        return _normalize_tick_reflections(chosen)
    except Exception:
        return []


def enrich_phase_results(result):
    """Attach consistent display labels and reflection labels after fitting."""
    tt = result.get("tt", []) or []
    tt_min = min(tt) if tt else 5.0
    tt_max = max(tt) if tt else 90.0
    wavelength = result.get("wavelength", 1.54056)
    phases = result.get("phase_results", []) or []
    for index, phase in enumerate(phases):
        phase["display_label"] = phase_display_label(phase, index=index)
        phase["legend_label"] = phase_legend_label(
            phase, index=index, show_weight=len(phases) > 1)
        phase["tick_label"] = phase_tick_label(phase, index=index)
        phase["tick_reflections"] = reflection_labels_for_phase(
            phase, wavelength, tt_min, tt_max)
    return result


def export_file_prefix(metadata):
    """Build a filesystem-safe ``YYYY-MM-DD_sample`` export prefix."""
    raw_date = metadata.get("analysis_date") or metadata.get("output_date")
    if isinstance(raw_date, (datetime, date)):
        date_text = raw_date.strftime("%Y-%m-%d")
    else:
        date_text = str(raw_date or "").strip()
        match = re.search(r"(\d{4})[-/]?(\d{2})[-/]?(\d{2})", date_text)
        date_text = "-".join(match.groups()) if match else date.today().isoformat()

    sample = str(metadata.get("sample_id") or "Sample").strip()
    sample = os.path.splitext(sample)[0]
    sample = unicodedata.normalize("NFKC", sample)
    sample = re.sub(r"[^\w.-]+", "_", sample, flags=re.UNICODE)
    sample = re.sub(r"[_-]{2,}", "_", sample).strip("_.-") or "Sample"
    return f"{date_text}_{sample}"
