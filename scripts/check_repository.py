"""Standard-library checks for repository navigation; never launches a solver."""
from pathlib import Path
import argparse
import json
import re
import tomllib
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
LIFECYCLE = {'planned', 'in_progress', 'complete', 'closed'}
VALIDATION = {'not_assessed', 'not_applicable', 'passed', 'failed', 'conditional', 'scenario_only', 'mixed'}


def check_repository(root=ROOT):
    root = Path(root).resolve()
    errors, examples, papers = [], {}, []
    for path in sorted((root/'examples').rglob('example.toml')):
        try:
            item = tomllib.loads(path.read_text())
            ident = item['id']
            if ident in examples:
                errors.append(f'Duplicate example ID: {ident}')
            examples[ident] = item
            if path.parent != root/'examples'/ident or item['area'] != ident.split('/')[0]:
                errors.append(f'Example ID/area differs from its folder: {path.relative_to(root)}')
            if item['schema_version'] != 1 or item['status'] not in LIFECYCLE:
                errors.append(f'Unsupported metadata version/status: {ident}')
            for key in ('readme', 'report'):
                if not (path.parent/item[key]).is_file():
                    errors.append(f'Missing {key}: {ident}')
            for key in ('reference_files', 'protocol_files'):
                for name in item[key]:
                    if not (path.parent/name).is_file():
                        errors.append(f'Missing {key} entry: {ident}/{name}')
            for key in ('run_module', 'check_module'):
                module = item[key]
                if module and not (root/Path(*module.split('.'))).with_suffix('.py').is_file():
                    errors.append(f'Missing {key}: {module}')
            for key, state in item['validation'].items():
                if key not in ('scope', 'note') and state not in VALIDATION:
                    errors.append(f'Invalid validation state: {ident}/{key}')
        except (KeyError, ValueError, TypeError) as error:
            errors.append(f'Invalid example metadata {path.relative_to(root)}: {error}')

    visited, active = set(), set()
    def visit(ident):
        if ident in active:
            errors.append(f'Cyclic example prerequisite: {ident}')
            return
        if ident in visited:
            return
        active.add(ident)
        for other in examples[ident].get('prerequisites', []):
            if other not in examples:
                errors.append(f'Unknown prerequisite {other} in {ident}')
            else:
                visit(other)
        active.remove(ident)
        visited.add(ident)
    for ident in examples:
        visit(ident)

    for path in sorted((root/'literature/sources').glob('*/metadata.json')):
        item = json.loads(path.read_text())
        papers.append(item['id'])
        if item['id'] != path.parent.name:
            errors.append(f'Paper ID differs from folder: {path.relative_to(root)}')
        for name in item['files']:
            if not (path.parent/name).is_file():
                errors.append(f'Missing paper evidence: {item["id"]}/{name}')

    documents = [root/'README.md', root/'ROADMAP.md']
    for folder in ('docs', 'examples', 'literature', 'templates', 'tests'):
        documents.extend((root/folder).rglob('*.md'))
    checked = 0
    for path in documents:
        relative = path.relative_to(root)
        if 'inbox' in relative.parts or (relative.parts[:2] == ('docs', 'history') and path.name != 'README.md'):
            continue  # Verbatim historical documents retain their original links.
        checked += 1
        for target in re.findall(r'\]\(([^\s)]+)\)', path.read_text()):
            if re.match(r'^[A-Za-z][A-Za-z+.-]*:', target) or target.startswith('#'):
                continue
            resolved = (path.parent/unquote(target.split('#')[0].strip('<>'))).resolve()
            if not resolved.is_relative_to(root):
                errors.append(f'File link leaves checkout: {relative}: {target}')
                continue
            if resolved.relative_to(root).parts[0] in ('results', '.local'):
                continue  # Runtime evidence may be unavailable in a source-only clone.
            if not resolved.exists():
                errors.append(f'Broken file link: {relative}: {target}')

    return {'examples': len(examples), 'paper_folders': len(papers),
            'markdown_files_checked': checked, 'errors': errors}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optional JSON record; no scientific files are updated.')
    args = parser.parse_args()
    report = check_repository()
    text = json.dumps(report, indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end='')
    raise SystemExit(bool(report['errors']))


if __name__ == '__main__':
    main()
