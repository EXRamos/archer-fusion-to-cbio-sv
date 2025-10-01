# Archer → cBioPortal SV — Prefect Project

This scaffold wraps your existing `archer_fusions_to_cbio_sv_v2.py` converter
in a Prefect flow and deployment, following the structure used by
`CBIIT/ChildhoodCancerDataInitiative-cBioPortal-Workflows` (folders: `src/`, `workflows/`, `prefect.yaml`, etc.).

## Project layout

```text
archer_prefect_project/
├── README.md
├── prefect.yaml
├── requirements.txt
├── src/
│   └── archer_cbio/
│       ├── __init__.py
│       ├── runner.py
│       └── bin/
│           └── archer_fusions_to_cbio_sv_v2.py
└── workflows/
    └── archer_to_cbio_flow.py
```

## Step-by-step: make your script Prefect-ready

1. **Put your code in `src/`**  
   We moved your converter into `src/archer_cbio/bin/archer_fusions_to_cbio_sv_v2.py` so it can be run as a module.

2. **Add a lightweight Python wrapper**  
   `src/archer_cbio/runner.py` exposes `run_archer_to_cbio(script_args)` which invokes the script via `python -m archer_cbio.bin.archer_fusions_to_cbio_sv_v2 ...`.

3. **Create a Prefect flow**  
   `workflows/archer_to_cbio_flow.py` defines `@flow archer_to_cbio_flow(...)` with two tasks:
   - `ensure_output_dir()` to create the output path
   - `run_converter()` to call the runner

4. **Declare a deployment in `prefect.yaml`**  
   The `prefect.yaml` describes how to pull/install code on workers and registers a deployment that points to the flow's entrypoint.

5. **Install and test locally**
   ```bash
   cd archer_prefect_project
   python -m venv .venv && source .venv/bin/activate  # or .venv\Scripts\activate on Windows
   pip install -r requirements.txt
   # quick local test (adjust args to match your script's CLI):
   python workflows/archer_to_cbio_flow.py
   ```

6. **Start a worker (in another terminal)**
   ```bash
   prefect worker start --pool default-process
   ```

7. **Deploy the flow**
   ```bash
   prefect deploy --name archer-fusions-to-cbio-sv
   # or specify by name if multiple deployments exist in prefect.yaml
   ```

8. **Run from UI or CLI**
   - UI: Navigate to your deployment and click **Run**.
   - CLI: `prefect deployment run 'archer-fusions-to-cbio-sv' --param script_args="['--input','/path/*.txt','--output','./out/data_sv.txt']"`

## Parameters

- `output_dir` – folder to create before running the converter.
- `script_args` – free-form list of CLI args passed to the original converter (keeps you flexible as the script evolves).
