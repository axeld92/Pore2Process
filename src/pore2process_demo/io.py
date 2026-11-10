"""Shared strict JSON and file hashing; no simulator initialization."""
import hashlib
import json
from pathlib import Path
from .paths import project_path

def read_json(path):
    return json.loads(project_path(Path(path).absolute()).read_text())

def write_json(path, value):
    path = project_path(Path(path).absolute())
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha(path):
    return hashlib.sha256(project_path(Path(path).absolute()).read_bytes()).hexdigest()
