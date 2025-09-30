from prefect import flow, task, get_run_logger
from pathlib import Path
from typing import Sequence, Optional, List
from archer_cbio.runner import run_archer_to_cbio

@task
def ensure_output_dir(path: str) -> str:
    p = Path(path).expanduser().resolve()
    p.mkdir(parents=True, exist_ok=True)
    return str(p)

@task
def run_converter(script_args: Sequence[str]) -> int:
    return run_archer_to_cbio(script_args)

@flow(name="archer-fusions-to-cbio-sv")
def archer_to_cbio_flow(
    output_dir: str = "./out",
    script_args: Optional[List[str]] = None
) -> str:
    """
    Prefect flow to execute the Archer→cBioPortal SV converter.

    Parameters
    ----------
    output_dir : str
        Directory to create/ensure exists. Useful when the converter expects an output path.
    script_args : list[str] | None
        Raw CLI arguments to pass directly to the converter script. This allows
        full control even if the underlying script's interface changes.
    """
    logger = get_run_logger()
    out_dir = ensure_output_dir(output_dir)
    logger.info(f"Output directory: {out_dir}")
    args = script_args or []
    rc = run_converter(args)
    if rc != 0:
        raise RuntimeError(f"Converter exited with non-zero status: {rc}")
    return out_dir

if __name__ == "__main__":
    # Example local run (adjust arguments to match your script's CLI):
    # archer_to_cbio_flow(output_dir="./out", script_args=["--input", "tests/data/*.txt", "--output", "./out/data_sv.txt'])
    archer_to_cbio_flow()
