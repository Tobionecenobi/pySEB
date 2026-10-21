"""Reuse expensive F127 scene exports across notebook-only changes."""
from functools import lru_cache
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import sys

from docs._prepare_cache import source_inputs


@lru_cache(None)
def implementation_key(root):
    digest = hashlib.sha256(sys.version.encode())
    for package in ('numpy', 'scipy', 'sympy', 'plotly'):
        digest.update(f'{package}={version(package)}'.encode())
    for name, content in source_inputs(root):
        if Path(name).name == 'compare-f127-coronas.qmd':
            # Track the export cell itself, without coupling it to the equation.
            content = json.dumps([cell for cell in json.loads(content)
                                  if 'scene_key(ROOT' in cell]).encode()
        elif name.endswith('.qmd') and Path(name).name not in (
                '_f127-model.qmd', '_f127-view.qmd'):
            continue
        digest.update(name.encode())
        digest.update(hashlib.sha256(content).digest())
    import pyseb
    native = Path(pyseb._pyseb.__file__)
    digest.update(hashlib.sha256(native.read_bytes()).digest())
    return digest.hexdigest()


def scene_key(root, model, arms, parameters, options, title):
    # Include all public option values, including defaults.
    settings = {name: getattr(options, name) for name in dir(options)
                if not name.startswith('_') and not callable(getattr(options, name))}
    payload = [implementation_key(root), model, arms, parameters, settings, title]
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def output_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reuse_scene(path, key):
    manifest = path.with_suffix('.cache.json')
    try:
        saved = json.loads(manifest.read_text())
        outputs = (path, path.with_suffix('.html'), path.parent / 'plotly.min.js')
        return saved['key'] == key and all(
            saved['outputs'][output.name] == output_hash(output) for output in outputs)
    except (OSError, ValueError, KeyError):
        return False


def record_scene(path, key):
    outputs = (path, path.with_suffix('.html'), path.parent / 'plotly.min.js')
    payload = {'key': key, 'outputs': {output.name: output_hash(output) for output in outputs}}
    manifest = path.with_suffix('.cache.json')
    temporary = manifest.with_suffix('.tmp')
    temporary.write_text(json.dumps(payload, indent=2) + '\n')
    temporary.replace(manifest)
