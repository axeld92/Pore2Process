"""Exercise native support-driver validation without substituting a Python solve."""
import subprocess
from datetime import datetime, timezone
from uuid import uuid4
from coupling.run import sha
from hydraulic.build import ROOT
from hydraulic.run import write_json
from filtration.build import BINARY


def main():
    directory=ROOT / "results/filtration/input_checks" / uuid4().hex[:8]
    directory.mkdir(parents=True)
    cases=(("nonfinite_radius",["--radius-m","nan"]),
           ("odd_grid",["--resolution","17"]),
           ("zero_pressure",["--pressure-drop-Pa","0"]),
           ("radius_outside_case",["--radius-m","6e-6"]),
           ("invalid_domain",["--length-m","1"]),
           ("trailing_numeric_text",["--radius-m","5e-6oops"]),
           ("unconverged",["--max-steps","1"]))
    results=[]
    for name,options in cases:
        work=directory / name;work.mkdir()
        command=[str(BINARY),*options]
        result=subprocess.run(command,cwd=work,text=True,capture_output=True,timeout=30)
        (work / "stdout.log").write_text(result.stdout+result.stderr)
        if result.returncode!=1 or "Filtration support reference failed:" not in result.stderr:
            raise RuntimeError(f"Native invalid/unconverged input did not fail cleanly: {name}")
        results.append({"name":name,"command":command,"returncode":result.returncode,
                        "error":result.stderr.strip()})
    report={"status":"passed","timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "binary_sha256":sha(BINARY),"cases":results}
    write_json(directory / "summary.json",report)
    print(f"Passed {len(results)} native invalid/unconverged cases: {directory / 'summary.json'}")


if __name__=="__main__":
    main()
