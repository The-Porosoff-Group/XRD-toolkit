"""Preserve GSAS-II's native fitting state without a parameter allowlist.

The human-readable snapshot is deliberately separate from display results: it
retains native units, full precision, refinement flags and project controls.
The accompanying GPX is the authoritative, reloadable native representation.
"""

import base64
import hashlib
import math
import platform
from collections.abc import Mapping
from importlib import metadata
from pathlib import Path

import numpy as np


def native_value(value, notes=None, path='$', _active=None):
    """Copy a native tree to strict JSON, recording unusual representations.

    NumPy arrays become lists; non-finite numbers and non-string mapping keys
    use explicit tagged representations. Unknown GSAS objects retain their
    class, text and attributes rather than silently becoming null or vanishing.
    """
    notes = notes if notes is not None else []
    active = _active if _active is not None else set()
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if math.isfinite(value):
            return value
        return {'__native_type__': 'float', 'value': str(value)}
    if isinstance(value, np.generic):
        return native_value(value.item(), notes, path, active)
    identity = id(value)
    if identity in active:
        notes.append(f'{path}: circular reference retained as a marker; use GPX.')
        return {'__native_type__': 'circular_reference'}
    active.add(identity)
    try:
        if isinstance(value, np.ma.MaskedArray):
            return {
                '__native_type__': 'masked_array',
                'data': native_value(value.data, notes, path + '.data', active),
                'mask': native_value(np.ma.getmaskarray(value), notes,
                                     path + '.mask', active),
                'fill_value': native_value(value.fill_value, notes,
                                           path + '.fill_value', active),
            }
        if isinstance(value, np.ndarray):
            return native_value(value.tolist(), notes, path, active)
        if isinstance(value, Mapping):
            if all(isinstance(key, str) for key in value):
                return {
                    key: native_value(item, notes, f'{path}/{key}', active)
                    for key, item in value.items()
                }
            return {
                '__native_type__': 'mapping',
                'entries': [
                    {'key': native_value(key, notes, f'{path}/key/{i}', active),
                     'value': native_value(item, notes, f'{path}/value/{i}', active)}
                    for i, (key, item) in enumerate(value.items())
                ],
            }
        if isinstance(value, (list, tuple)):
            return [native_value(item, notes, f'{path}/{i}', active)
                    for i, item in enumerate(value)]
        if isinstance(value, (set, frozenset)):
            return {'__native_type__': type(value).__name__,
                    'values': native_value(sorted(value, key=repr), notes,
                                           path + '/values', active)}
        if isinstance(value, (bytes, bytearray)):
            return {'__native_type__': 'bytes',
                    'base64': base64.b64encode(value).decode('ascii')}
        kind = f'{type(value).__module__}.{type(value).__qualname__}'
        notes.append(f'{path}: {kind} preserved as text and attributes; '
                     'the GPX retains its native type.')
        return {
            '__native_type__': kind,
            'text': str(value),
            'attributes': native_value(getattr(value, '__dict__', {}), notes,
                                       path + '/attributes', active),
        }
    finally:
        active.remove(identity)


def project_snapshot(project, notes=None, compact=False):
    """Capture native project data, including future/unrecognized parameters.

    Stage snapshots omit only repeated scan/reflection arrays and covariance;
    initial/final snapshots include those fields in full. Every omission is
    named so a workbook reader can distinguish it from an absent parameter.
    """
    data = project.data
    excluded = []
    if compact:
        data = dict(data)
        if 'Covariance' in data:
            data.pop('Covariance')
            excluded.append('data/Covariance (retained in full final_project)')
        for name, content in list(data.items()):
            if str(name).startswith('PWDR ') and isinstance(content, Mapping):
                content = dict(content)
                for key in ('data', 'Reflection Lists'):
                    if key in content:
                        content.pop(key)
                        excluded.append(f'data/{name}/{key} '
                                        '(retained in initial/final_project)')
                data[name] = content
    return {
        'data': native_value(data, notes),
        'names': native_value(getattr(project, 'names', []), notes),
        'excluded_paths': excluded,
    }


def software_provenance(gsas_module, backend_file):
    """Identify the installed code without requiring network/git access."""
    result = {'python': platform.python_version(), 'numpy': np.__version__,
              'platform': platform.platform()}
    for package in ('GSAS-II', 'GSASII', 'gsas2pkg', 'scipy', 'pymatgen'):
        try:
            result[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            pass
    for label, filename in (
            ('gsas_scriptable', getattr(gsas_module, '__file__', None)),
            ('toolkit_backend', backend_file),
            ('native_capture', __file__)):
        if filename:
            path = Path(filename)
            result[label] = {'filename': path.name}
            try:
                result[label]['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            except OSError as exc:
                result[label]['read_error'] = str(exc)
    for name in ('__version__', '__revision__', 'git_version'):
        value = getattr(gsas_module, name, None)
        if value is not None and not callable(value):
            result[f'gsas_{name.strip("_")}'] = str(value)
    # Installed GSAS-II distributions need not expose package metadata or
    # __version__. Their bundled version record identifies the source commit.
    gsas_path = getattr(gsas_module, 'GSASIIpath', None)
    get_saved_version = getattr(gsas_path, 'getSavedVersionInfo', None)
    if callable(get_saved_version):
        try:
            saved_version = get_saved_version()
            if saved_version is not None:
                result['gsas_version'] = {
                    name: native_value(getattr(saved_version, name))
                    for name in ('git_version', 'git_versiontag', 'git_tags',
                                 'git_prevtaggedversion', 'git_prevtags')
                    if hasattr(saved_version, name)
                }
        except Exception as exc:
            result['gsas_version_read_error'] = str(exc)
    return result


class NativeFitRecorder:
    """Keep snapshots independent of later phase-isolation/display mutations."""

    def __init__(self, submitted_settings, gsas_module, backend_file):
        self.notes = []
        self.payload = {
            'schema_version': 1,
            'capture_notes': self.notes,
            'value_convention': ('Native GSAS-II units and full precision. '
                                 'List positions and refinement flags are retained. '
                                 'NumPy arrays are represented as JSON lists.'),
            'submitted_settings': native_value(submitted_settings, self.notes),
            'effective_settings': {},
            'input_files': {},
            'software': software_provenance(gsas_module, backend_file),
            'refinement_stages': [],
        }
        self.project_base64 = None

    def begin_stage(self, project, name, number, recipe):
        if 'initial_project' not in self.payload:
            self.payload['initial_project'] = project_snapshot(project, self.notes)
        stage = {
            'name': name, 'number': number,
            'recipe': native_value(recipe, self.notes),
            'before': project_snapshot(project, self.notes, compact=True),
        }
        self.payload['refinement_stages'].append(stage)
        return stage

    def finish(self, project, project_path):
        # This must happen before the backend clears cell flags or evaluates
        # one phase at a time: both actions alter the native fitted state.
        self.payload['final_project'] = project_snapshot(project, self.notes)
        try:
            project.save()
            contents = Path(project_path).read_bytes()
            self.project_base64 = base64.b64encode(contents).decode('ascii')
            self.payload['native_project'] = {
                'format': 'GSAS-II GPX',
                'sha256': hashlib.sha256(contents).hexdigest(),
                'bytes': len(contents),
                'capture_point': 'After fitting, before phase-isolation/display changes',
            }
        except Exception as exc:
            self.notes.append(f'Native GPX companion could not be saved: {exc}')
