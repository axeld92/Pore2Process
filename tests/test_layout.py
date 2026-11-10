"""Location-independent discovery and guarded historical source replay."""
import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from pore2process_demo.paths import ROOT, project_root
from pore2process_demo.provenance import current_source_path, implementation_files


def test_root_and_installed_launcher_work_outside_checkout(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert project_root() == ROOT
    python = ROOT / '.venv/bin/python'
    result = subprocess.run([str(python), '-m', 'examples.hydraulics.capillary_network.run', '--help'],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert '--output' in result.stdout


def test_old_source_names_resolve_without_old_top_level_packages():
    assert current_source_path('hydraulic/run.py') == ROOT/'examples/hydraulics/capillary_benchmark/run.py'
    assert current_source_path('compatibility/dwsim_runner/Program.cs') == ROOT/'src/pore2process_demo/runtime/dwsim_runner/Program.cs'
    assert not (ROOT/'hydraulic').exists()
    with pytest.raises(ValueError):current_source_path('../outside.py')


def test_new_snapshots_include_extracted_shared_helpers():
    files = implementation_files([ROOT/'examples/filtration/operating_map/run.py'])
    assert ROOT/'src/pore2process_demo/io.py' in files
    assert ROOT/'src/pore2process_demo/runtime/openlb.py' in files
    assert ROOT/'src/pore2process_demo/runtime/hydraulics.py' in files
    assert ROOT/'docs/history/reorganization.json' in files


@pytest.fixture
def isolated_provenance(tmp_path, monkeypatch):
    from pore2process_demo import provenance as p
    digest = lambda b: hashlib.sha256(b).hexdigest()
    current=tmp_path/'new.py';current.write_bytes(b'current relocated source')
    old=b'frozen historical source';expected=digest(old)
    archived=tmp_path/'.local/reorganization/source_cache/revision/old.py'
    archived.parent.mkdir(parents=True);archived.write_bytes(old)
    manifest=tmp_path/'migration.json'
    manifest.write_text(json.dumps({'paths':{'old.py':'new.py'},'before_revision':'revision',
        'sources':{'old.py':{'before_sha256':expected,'current_guards':{'new.py':digest(current.read_bytes())}}}}))
    monkeypatch.setattr(p,'ROOT',tmp_path);monkeypatch.setattr(p,'MANIFEST',manifest);p._manifest.cache_clear()
    try:yield p,current,archived,expected
    finally:p._manifest.cache_clear()


def test_historical_hash_requires_both_archived_and_current_source(isolated_provenance):
    p,current,archived,expected=isolated_provenance
    assert p.verified_source_path('old.py',expected) == archived
    current.write_bytes(b'changed scientific implementation')
    with pytest.raises(ValueError,match='Changed relocated implementation'):
        p.verified_source_path('old.py',expected)


def test_corrupted_historical_bytes_fail(isolated_provenance):
    p,current,archived,expected=isolated_provenance
    archived.write_bytes(b'corrupted archive')
    with pytest.raises(ValueError,match='Changed archived source'):
        p.verified_source_path('old.py',expected)
