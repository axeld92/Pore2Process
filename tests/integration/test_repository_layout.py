"""Resource moves must retain provenance and fail on unreviewed edits."""
import hashlib
import json
from pathlib import Path
import runpy

import pytest

from pore2process_demo.paths import ROOT, project_path
from pore2process_demo.provenance import implementation_files, verified_source_path


def test_frozen_resource_aliases_have_exact_original_bytes():
    registry = json.loads((ROOT/'docs/history/resource_organization.json').read_text())
    for name in ('docs/film_particle_reference.json', 'literature/zhu2021_fig6_pdf_points.json',
                 'literature/alhamed2009_fig7_pixels.csv', 'examples/filtration/graded_support/protocol.json'):
        record = registry['sources'][name]
        actual = project_path(name)
        assert actual.is_file()
        assert hashlib.sha256(actual.read_bytes()).hexdigest() == record['before_sha256']
        assert verified_source_path(name, record['before_sha256']) == actual


def test_new_snapshots_include_case_metadata_and_both_registries():
    files = implementation_files([ROOT/'examples/filtration/operating_map/run.py'])
    assert ROOT/'examples/filtration/operating_map/example.toml' in files
    assert ROOT/'docs/history/reorganization.json' in files
    assert ROOT/'docs/history/resource_organization.json' in files


def test_later_reviewed_transition_preserves_old_source_and_rejects_edits(tmp_path, monkeypatch):
    from pore2process_demo import provenance as p
    digest = lambda value: hashlib.sha256(value).hexdigest()
    original, middle, current = b'first source', b'after first move', b'after second move'
    final = tmp_path/'new.py';final.write_bytes(current)
    history = tmp_path/'docs/history';history.mkdir(parents=True)
    first = history/'reorganization.json'
    first.write_text(json.dumps({'before_revision':'first','paths':{'old.py':'middle.py'},
        'sources':{'old.py':{'before_sha256':digest(original),'current_guards':{'middle.py':digest(middle)}}}}))
    (history/'resource_organization.json').write_text(json.dumps({'before_revision':'second',
        'paths':{'middle.py':'new.py'},'sources':{'middle.py':{'before_sha256':digest(middle),
        'current_guards':{'new.py':digest(current)}}}}))
    archive=tmp_path/'.local/reorganization/source_cache/first/old.py'
    archive.parent.mkdir(parents=True);archive.write_bytes(original)
    monkeypatch.setattr(p,'ROOT',tmp_path);monkeypatch.setattr(p,'MANIFEST',first);p._manifest.cache_clear()
    try:
        assert p.verified_source_path('old.py',digest(original)) == archive
        final.write_bytes(b'unreviewed implementation')
        with pytest.raises(ValueError,match='Changed relocated implementation'):
            p.verified_source_path('old.py',digest(original))
    finally:
        p._manifest.cache_clear()


def test_catalogue_and_prerequisites_are_consistent_without_running_examples():
    check = runpy.run_path(str(ROOT/'scripts/check_repository.py'))['check_repository']
    report = check(ROOT)
    assert report['examples'] >= 16 and report['paper_folders'] >= 13
    assert report['errors'] == []
