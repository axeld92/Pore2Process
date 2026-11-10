"""Install the official pinned .NET 10 SDK privately, with publisher SHA512."""
from pathlib import Path
import hashlib
import json
import platform
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / ".local/compatibility"
METADATA_URL = "https://builds.dotnet.microsoft.com/dotnet/release-metadata/10.0/releases.json"
VERSION = "10.0.403"
SHA512 = "5893c66b48ea3674c7da18501059d606e96f3672fabe4d5e41762e692c7c83df3ee0b2b9dc5a0786191da2dc413d7335158101c14c935ab217f31ac2529770f3"


def main():
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        raise SystemExit("This pinned SDK is for Apple Silicon macOS; select the official SDK for this machine before adapting it.")
    (CACHE / "downloads").mkdir(parents=True, exist_ok=True)
    metadata_path = CACHE / "dotnet-releases.json"
    if not metadata_path.exists():
        with urllib.request.urlopen(METADATA_URL,timeout=30) as response:
            metadata_path.write_bytes(response.read())
    metadata = json.loads(metadata_path.read_text())
    sdk = next(sdk for release in metadata["releases"] for sdk in release.get("sdks",[release.get("sdk",{})]) if sdk.get("version") == VERSION)
    asset = next(file for file in sdk["files"] if file["rid"] == "osx-arm64" and file["name"].endswith(".tar.gz"))
    if asset["hash"] != SHA512:
        raise RuntimeError("SDK publisher metadata differs from the tested pinned SHA512.")
    archive = CACHE / "downloads" / Path(asset["url"]).name
    if not archive.exists():
        print(f"Downloading official .NET SDK {VERSION} for Apple Silicon...",flush=True)
        with urllib.request.urlopen(asset["url"],timeout=45) as response, archive.with_suffix(".part").open("wb") as output:
            while chunk := response.read(1024*1024):
                output.write(chunk)
        archive.with_suffix(".part").replace(archive)
    with archive.open("rb") as source_file:
        digest = hashlib.file_digest(source_file,"sha512").hexdigest()
    if digest != SHA512:
        raise RuntimeError("Publisher SHA512 mismatch; refusing to extract .NET.")
    destination = CACHE / "dotnet"
    if not destination.exists():
        with tempfile.TemporaryDirectory(prefix="extract-", dir=CACHE) as staging:
            payload = Path(staging) / "payload"
            payload.mkdir()
            with tarfile.open(archive) as package:
                package.extractall(payload,filter="data")
            payload.rename(destination)
    if not (destination / "dotnet").is_file():
        raise RuntimeError(f"Incomplete SDK cache at {destination}; inspect it before removing/recreating it.")
    record = {"sdk_version":VERSION,"metadata_url":METADATA_URL,**asset,"sha512":digest,"directory":str(destination)}
    (CACHE / "dotnet_distribution.json").write_text(json.dumps(record,indent=2)+"\n")
    print(f"Private SDK ready: {destination / 'dotnet'}",flush=True)


if __name__ == "__main__":
    main()
