import subprocess
import sys
from pathlib import Path
from typing import Sequence, Optional

def run_archer_to_cbio(script_args: Optional[Sequence[str]] = None, python_executable: str = sys.executable) -> int:
    """Run the Archer → cBioPortal SV converter as a Python module via subprocess.

    Parameters
    ----------
    script_args : Sequence[str], optional
        Raw CLI arguments to pass through to the original script.
        Example: ["--input", "input_dir/*.txt", "--output", "out/data_sv.txt"]
    python_executable : str
        Path to the Python interpreter to use. Defaults to current interpreter.

    Returns
    -------
    int
        The process return code (0 indicates success).
    """
    # We run the script using -m so it can still reside inside the package.
    module = "archer_cbio.bin.archer_fusions_to_cbio_sv_v2"
    cmd = [python_executable, "-m", module]
    if script_args:
        cmd.extend(script_args)
    print("[runner] Executing:", " ".join(cmd))
    proc = subprocess.run(cmd, check=False)
    return proc.returncode
