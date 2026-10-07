"""Helpers for producing standards-compliant JSON values."""

import math

import numpy as np


def json_safe_value(value):
    """Recursively replace non-finite numbers with JSON ``null`` values."""
    if isinstance(value, dict):
        return {key: json_safe_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe_value(value.tolist())
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return number if math.isfinite(number) else None
    return value
