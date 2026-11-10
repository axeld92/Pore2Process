"""Shared strict JSON and file hashing; no simulator initialization."""
import hashlib
import json
from pathlib import Path

def read_json(path):
    return json.loads(Path(path).read_text())

def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
