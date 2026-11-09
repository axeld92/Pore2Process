"""Download official, pinned tool distributions into this project's ignored cache."""
from pathlib import Path
import hashlib
import json
import platform
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / ".local/compatibility"
SOURCES = {
    "openlb": {
        "version": "1.9r0",
        "url": "https://zenodo.org/records/17899765/files/olb-1.9r0.tgz?download=1",
        "filename": "olb-1.9r0.tgz",
        "sha256": "e4d2421b50643482036917c4416145f7d6b495492f3afac64e048212361c8b3e",
    },
    "dwsim": {
        "version": "10.2.15",
        "url": "https://github.com/DanWBR/dwsim10/releases/download/v10.2.15/dwsim-10.2.15-osx-arm64.tar.gz",
        "filename": "dwsim-10.2.15-osx-arm64.tar.gz",
        "sha256": "614a214178dc2c1c46156f355b2b8aeca30d66758ca62be11b7cc67f9a95fde6",
    },
}


def main():
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        raise SystemExit("This pinned compatibility setup targets Apple Silicon macOS; choose official distributions for this machine before adapting it.")
    downloads = CACHE / "downloads"
    downloads.mkdir(parents=True,exist_ok=True)
    records = {}
    for name,source in SOURCES.items():
        archive = downloads / source["filename"]
        if not archive.exists():
            print(f"Downloading official {name} {source['version']}...",flush=True)
            request = urllib.request.Request(source["url"],headers={"User-Agent":"Pore2Process-compatibility/0.1"})
            partial = archive.with_suffix(".part")
            with urllib.request.urlopen(request,timeout=45) as response, partial.open("wb") as output:
                while chunk := response.read(1024*1024):
                    output.write(chunk)
            partial.replace(archive)
        with archive.open("rb") as source_file:
            digest = hashlib.file_digest(source_file,"sha256").hexdigest()
        if source["sha256"] and digest != source["sha256"]:
            raise RuntimeError(f"{name} SHA256 mismatch; refusing to extract.")
        destination = CACHE / name
        if not destination.exists():
            with tempfile.TemporaryDirectory(prefix="extract-", dir=CACHE) as staging:
                payload = Path(staging) / "payload"
                payload.mkdir()
                with tarfile.open(archive) as package:
                    package.extractall(payload,filter="data")
                payload.rename(destination)
        expected = destination / ("olb-1.9r0/Makefile" if name == "openlb" else "DWSIM.Automation.FluentAPI.dll")
        if not expected.is_file():
            raise RuntimeError(f"Incomplete distribution cache at {destination}; inspect it before removing/recreating it.")
        records[name] = {**source,"sha256":digest,"bytes":archive.stat().st_size,"directory":str(destination)}
        print(f"{name}: {archive.stat().st_size} bytes; SHA256 {digest}; {destination}",flush=True)
    (CACHE / "distributions.json").write_text(json.dumps(records,indent=2)+"\n")


if __name__ == "__main__":
    main()
