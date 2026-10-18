"""Invalidate Quarto notebook caches when their external inputs change."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = 1
PYTHON_CELLS = re.compile(r'^```(?:\{python\}|python)[ \t]*\n(.*?)^```', re.M | re.S)
SOURCE_SUFFIXES = {'.py', '.cpp', '.c', '.h', '.hpp', '.yaml', '.toml', '.cmake', '.in'}


def source_inputs(root):
    """Include calculation code and data, but exclude page prose and layout."""
    files = {root / 'CMakeLists.txt', root / 'pyproject.toml'}
    workflow = root / '.github/workflows/documentation.yml'
    if workflow.is_file():
        yield workflow.relative_to(root).as_posix(), workflow.read_bytes()
    for directory in ('seb', 'seb-symbolic', 'pyseb', 'models', 'third_party', 'cmake'):
        files.update(path for path in (root / directory).rglob('*')
                     if path.is_file() and (path.suffix in SOURCE_SUFFIXES
                                           or path.name == 'CMakeLists.txt'))
    docs = root / 'docs'
    files.update(docs.rglob('*.py'))
    files.update((docs / 'tutorials/cpp').rglob('*'))
    files.update((docs / 'tutorials/data').rglob('*'))
    files.update((docs / 'tutorials/models').rglob('*'))
    for path in sorted(files):
        if path.is_relative_to(docs) and any(part in ('_site', '_freeze', 'generated')
                                             for part in path.relative_to(docs).parts):
            continue
        if path.is_file() and not any(part.startswith('.') or part == '__pycache__'
                                      for part in path.relative_to(root).parts):
            yield path.relative_to(root).as_posix(), path.read_bytes()
    for path in sorted(docs.rglob('*.qmd')):
        if any(part.startswith('.') or part in ('_site', '_freeze')
               for part in path.relative_to(docs).parts):
            continue
        # Keep cached notebook outputs and their generated downloads in the
        # same calculation context, including when reverting to older code.
        cells = PYTHON_CELLS.findall(path.read_text())
        if cells:
            yield path.relative_to(root).as_posix(), json.dumps(cells).encode()


def runtime_description():
    python = os.environ.get('QUARTO_PYTHON', sys.executable)
    environment = subprocess.check_output([
        python, '-c',
        'import platform,subprocess,sys; '
        'print(sys.version); print(platform.system(), platform.machine()); '
        'print(subprocess.check_output([sys.executable,"-m","pip","freeze"],text=True))',
    ], text=True)
    quarto = subprocess.check_output(['quarto', '--version'], text=True)
    native = []
    for command in (['gsl-config', '--version'], ['c++', '--version']):
        if shutil.which(command[0]):
            result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, check=False)
            native.append([command[0], result.returncode, result.stdout])
    return environment + '\nQuarto ' + quarto + json.dumps(native)


def fingerprint(root, runtime):
    digest = hashlib.sha256(f'docs-cache-{SCHEMA}\n{runtime}'.encode())
    for name, content in source_inputs(root):
        digest.update(name.encode() + b'\0')
        digest.update(hashlib.sha256(content).digest())
    return digest.hexdigest()


def metadata_path(root):
    return root / 'docs/.jupyter_cache/pyseb-inputs.json'


def write_metadata(root, data):
    path = metadata_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def prepare(root, key):
    path = metadata_path(root)
    try:
        previous = json.loads(path.read_text())
    except (OSError, ValueError):
        previous = {}
    missing = any(not (root / name).is_file()
                  for name in previous.get('generated_outputs', []))
    if previous.get('fingerprint') == key and not missing:
        print('Documentation: reuse notebook results (calculation inputs unchanged).')
        return False
    docs = root / 'docs'
    cache_directories = {docs / '.jupyter_cache'}
    cache_directories.update(path.parent / '.jupyter_cache' for path in docs.rglob('*.qmd')
                             if '_site' not in path.parts and '_freeze' not in path.parts)
    for directory in cache_directories:
        if directory.is_dir():
            shutil.rmtree(directory)
    write_metadata(root, {'fingerprint': key, 'generated_outputs': []})
    reason = 'generated downloads are missing' if missing else 'calculation inputs or environment changed'
    print(f'Documentation: refresh notebook results ({reason}).')
    return True


def record_outputs(root):
    path = metadata_path(root)
    if not path.is_file():
        return
    metadata = json.loads(path.read_text())
    directory = root / 'docs/tutorials/generated'
    metadata['generated_outputs'] = [path.relative_to(root).as_posix()
                                     for path in sorted(directory.rglob('*')) if path.is_file()]
    write_metadata(root, metadata)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--key', action='store_true', help='Print the CI cache identity without modifying files')
    parser.add_argument('--record-outputs', action='store_true')
    args = parser.parse_args()
    if args.record_outputs:
        record_outputs(ROOT)
    else:
        key = fingerprint(ROOT, runtime_description())
        if args.key:
            print(key)
        else:
            prepare(ROOT, key)


if __name__ == '__main__':
    main()
