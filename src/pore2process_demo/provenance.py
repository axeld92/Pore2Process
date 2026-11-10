"""Verify frozen pre-move sources without rewriting accepted run records.

Current code must match the reviewed relocation manifest before an archived
source version can satisfy a historical hash. Numerical/artifact checks remain
with each example's existing validator.
"""
from functools import lru_cache
import json
from pathlib import Path
import subprocess

from .paths import ROOT
from .io import sha

MANIFEST = ROOT / 'docs/history/reorganization.json'

@lru_cache(maxsize=1)
def _manifest():
    return json.loads(MANIFEST.read_text())


def current_source_path(name):
    path = Path(name)
    if path.is_absolute():
        try:
            path = path.relative_to(ROOT)
        except ValueError as error:
            raise ValueError('Source must belong to this project.') from error
    if '..' in path.parts:
        raise ValueError('Source path must not traverse outside the project.')
    relative = path.as_posix()
    mapping = _manifest()['paths']
    for old, new in sorted(mapping.items(), key=lambda item: -len(item[0])):
        if relative == old or relative.startswith(old + '/'):
            return ROOT / (new + relative[len(old):])
    return ROOT / path


def verified_source_path(name, expected_sha256):
    """Return bytes with the recorded hash, checking both sides of a relocation."""
    current = current_source_path(name)
    if current.is_file() and sha(current) == expected_sha256:
        return current
    record = _manifest()['sources'].get(str(name))
    if not record or record['before_sha256'] != expected_sha256:
        raise ValueError(f'Changed source without a verified relocation: {name}')
    for relative, digest in record['current_guards'].items():
        candidate = ROOT / relative
        if not candidate.is_file() or sha(candidate) != digest:
            raise ValueError(f'Changed relocated implementation: {relative}')
    revision = _manifest()['before_revision']
    target = ROOT / '.local/reorganization/source_cache' / revision / name
    if not target.exists():
        result = subprocess.run(['git', 'show', f'{revision}:{name}'], cwd=ROOT,
                                capture_output=True, timeout=15)
        if result.returncode:
            raise ValueError('Historical source revision is unavailable; fetch full Git history before replay.')
        import hashlib
        if hashlib.sha256(result.stdout).hexdigest() != expected_sha256:
            raise ValueError(f'Historical source digest differs: {name}')
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(result.stdout)
    if sha(target) != expected_sha256:
        raise ValueError(f'Changed archived source: {name}')
    return target


def implementation_files(paths):
    """Include shared helpers in all new source/input snapshots."""
    package = ROOT / 'src/pore2process_demo'
    files = [Path(p) if Path(p).is_absolute() else ROOT / p for p in paths]
    files += list(package.glob('*.py')) + list((package / 'runtime').glob('*.py'))
    files += [MANIFEST]
    return sorted({p for p in files if p.is_file()})


def source_manifest(paths):
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in implementation_files(paths)}
