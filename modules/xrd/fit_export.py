"""Lossless, inspectable fit provenance rows for XRD workbooks."""

import json
import os

import numpy as np


def json_value(value):
    """Copy fit inputs into plain JSON values without rounding numbers."""
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    return value


def fit_parameter_rows(result, metadata, method_label):
    """Flatten dictionaries; retain arrays as JSON and split oversized cells.

    Values are JSON text, preserving float precision, booleans, nulls and
    literal strings. Rejoin numbered parts before decoding a long value.
    Parameter paths use JSON Pointer escaping (~0 for ~ and ~1 for /).
    """
    payload = {
        'Export guide': {
            'format': 'Values are JSON text. Join Parts in order before decoding. '
                      'Parameter paths use JSON Pointer escaping.',
            'native_state': 'initial_project and final_project contain the native '
                            'GSAS-II project tree, including fixed/refined flags. '
                            'Final state precedes diagnostic phase isolation.',
            'reproduction': 'Use the companion GSAS-II .gpx project when available; '
                            'the workbook records inputs, settings and native state. '
                            'Use matching software versions for reproduction.',
            'method': method_label,
        },
        'Fit settings': result.get('fit_settings', {
            'availability': 'Not captured in this older result. Run the fit again '
                            'to record its inputs and settings.',
        }),
        'Figure settings at fit completion': metadata,
        'GSAS-II native parameters': result.get('gsas_native_parameters', {
            'availability': ('Native GSAS-II state was not captured. Run the fit '
                             'again with the updated toolkit.'
                             if 'gsas' in method_label.lower() else
                             'Not applicable: this fit did not use GSAS-II.'),
        }),
    }
    rows = []

    def visit(section, path, value):
        if isinstance(value, dict) and value:
            for key, item in value.items():
                escaped = str(key).replace('~', '~0').replace('/', '~1')
                visit(section, path + '/' + escaped, item)
            return
        encoded = json.dumps(json_value(value), ensure_ascii=True,
                             separators=(',', ':'))
        chunks = [encoded[i:i + 30000] for i in range(0, len(encoded), 30000)]
        for index, chunk in enumerate(chunks, 1):
            rows.append({
                'Section': section,
                'Parameter path': path or '/',
                'Type': ('null' if value is None else
                         'boolean' if isinstance(value, bool) else
                         'object' if isinstance(value, dict) else
                         'array' if isinstance(value, (list, tuple, np.ndarray)) else
                         'string' if isinstance(value, (str, os.PathLike)) else
                         'number'),
                'Part': index,
                'Parts': len(chunks),
                'Value (JSON)': chunk,
            })

    for section, value in payload.items():
        visit(section, '', value)
    return rows
