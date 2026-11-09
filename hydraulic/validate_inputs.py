"""Check native input validation and refusal to export an unconverged result."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]


def main():
    binary=ROOT / ".local/hydraulic/build/capillary3d"
    output=ROOT / "results/hydraulic/input_validation" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    output.mkdir(parents=True)
    records=[]
    cases=(("--resolution","12.5"),("--resolution","12extra"),("--resolution","7"),
           ("--resolution","9"),("--resolution","66"),("--pressure-drop-Pa","0"),
           ("--pressure-drop-Pa","nan"),("--pressure-drop-Pa","inf"),
           ("--tau","0.5"),("--steady-tolerance","0"),("--max-steps","0"),
           ("--unknown","1"),("--resolution",),("--resolution","8","--max-steps","1"),
           ("--radius-m","nan"),("--radius-m","inf"),("--radius-m","-1"),
           ("--radius-m","5e-6"),("--radius-m","0.0001"),("--radius-m","47e-6"),
           ("--viscosity-Pa-s","0"),("--viscosity-Pa-s","-1"),("--viscosity-Pa-s","nan"),
           ("--density-kg-m3","0"),("--density-kg-m3","nan"))
    for i,arguments in enumerate(cases):
        directory=output / str(i)
        directory.mkdir()
        result=subprocess.run([str(binary),*arguments],cwd=directory,text=True,capture_output=True,timeout=20)
        (directory / "stdout.log").write_text(result.stdout+result.stderr)
        if result.returncode != 1 or "Hydraulic benchmark failed:" not in result.stderr or (directory / "raw.json").exists():
            raise RuntimeError(f"Invalid/unconverged case did not fail cleanly: {arguments}")
        records.append({"arguments":arguments,"return_code":result.returncode,"exported_result":False,
                        "message":result.stderr.strip()})
    (output / "summary.json").write_text(json.dumps({"status":"passed","cases":records},indent=2)+"\n")
    print(f"Passed {len(records)} native invalid/unconverged checks: {output / 'summary.json'}")


if __name__ == "__main__":
    main()
