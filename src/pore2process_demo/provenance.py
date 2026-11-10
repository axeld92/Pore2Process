"""Replay source hashes through reviewed, ordered relocations.

Accepted records and source snapshots are never rewritten. An older source
guard may cross a later reviewed transition, never an unrelated source edit.
"""
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import subprocess

from .paths import ROOT

MANIFEST = ROOT / 'docs/history/reorganization.json'


@lru_cache(maxsize=1)
def _manifest():
    return json.loads(MANIFEST.read_text())


@lru_cache(maxsize=4)
def _additional_manifest(path):
    return json.loads(path.read_text())


def _manifests():
    records = [_manifest()]
    extra = ROOT / 'docs/history/resource_organization.json'
    if extra.is_file():
        records.append(_additional_manifest(extra))
    return records


def _relative(name):
    path = Path(name)
    if path.is_absolute():
        try:
            path = path.relative_to(ROOT)
        except ValueError as error:
            raise ValueError('Source must belong to this project.') from error
    if '..' in path.parts:
        raise ValueError('Source path must not traverse outside the project.')
    return path.as_posix()


def _move(name, manifest):
    for old, new in sorted(manifest['paths'].items(), key=lambda item: -len(item[0])):
        if name == old or name.startswith(old + '/'):
            return new + name[len(old):]
    return name


def _current(name, manifests):
    for manifest in manifests:
        name = _move(name, manifest)
    return ROOT / name


def current_source_path(name):
    return _current(_relative(name), _manifests())


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matching(name, digest, manifests, after=-1):
    for index, manifest in enumerate(manifests):
        record = manifest.get('sources', {}).get(name)
        if index > after and record and record['before_sha256'] == digest:
            return index, name, record
        name = _move(name, manifest)
    return None


def _verify_record(index, record, manifests):
    for relative, digest in record['current_guards'].items():
        path = _current(relative, manifests)
        if path.is_file() and _sha(path) == digest:
            continue
        later = _matching(relative, digest, manifests, after=index)
        if later is None:
            raise ValueError(f'Changed relocated implementation: {relative}')
        _verify_record(later[0], later[2], manifests)


def verified_source_path(name, expected_sha256):
    """Require exact recorded bytes and all applicable relocation guards."""
    name = _relative(name)
    manifests = _manifests()
    current = _current(name, manifests)
    if current.is_file() and _sha(current) == expected_sha256:
        return current
    match = _matching(name, expected_sha256, manifests)
    if match is None:
        raise ValueError(f'Changed source without a verified relocation: {name}')
    index, archived_name, record = match
    _verify_record(index, record, manifests)
    revision = record.get('before_revision', manifests[index]['before_revision'])
    target = ROOT / '.local/reorganization/source_cache' / revision / archived_name
    if not target.exists():
        result = subprocess.run(['git', 'show', f'{revision}:{archived_name}'], cwd=ROOT,
                                capture_output=True, timeout=15)
        if result.returncode:
            raise ValueError('Historical source revision is unavailable; fetch full Git history before replay.')
        if hashlib.sha256(result.stdout).hexdigest() != expected_sha256:
            raise ValueError(f'Historical source digest differs: {name}')
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(result.stdout)
    if _sha(target) != expected_sha256:
        raise ValueError(f'Changed archived source: {name}')
    return target


def implementation_files(paths):
    """Include resource resolvers and both registries in new run fingerprints."""
    package = ROOT / 'src/pore2process_demo'
    files = [current_source_path(p) for p in paths]
    for path in list(files):
        if path.is_relative_to(ROOT / 'examples'):
            for parent in path.parents:
                if parent == ROOT / 'examples':
                    break
                metadata = parent / 'example.toml'
                if metadata.is_file():
                    files.append(metadata)
                    break
    files += list(package.glob('*.py')) + list((package / 'runtime').glob('*.py'))
    files += [MANIFEST, ROOT / 'docs/history/resource_organization.json']
    return sorted({p for p in files if p.is_file()})


def source_manifest(paths):
    return {p.relative_to(ROOT).as_posix(): _sha(p) for p in implementation_files(paths)}
