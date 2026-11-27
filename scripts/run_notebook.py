"""Execute a notebook with the active Python kernel and save its outputs after each cell."""
import argparse
from tempfile import NamedTemporaryFile
from pathlib import Path
import nbformat
from nbclient import NotebookClient

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("notebook", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--preset", help="execute only the setup cells and this tagged experiment")
    args = parser.parse_args()
    path = args.notebook.resolve()
    out = args.output or path
    nb = nbformat.read(path, as_version=4)
    if args.preset:
        if out.resolve() == path:
            parser.error("--preset requires a separate --output notebook")
        if not any(args.preset in c.metadata.get("tags", []) for c in nb.cells):
            parser.error(f"unknown preset {args.preset}")
        nb.cells = [c for c in nb.cells if set(c.metadata.get("tags", [])) & {"setup", args.preset}]
    out.parent.mkdir(parents=True, exist_ok=True)
    for cell in nb.cells:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
            cell.metadata.pop("execution", None)
    nbformat.validate(nb)
    client = NotebookClient(nb, timeout=None, kernel_name="python3",
                            resources={"metadata": {"path": str(path.parent)}})
    def checkpoint():
        with NamedTemporaryFile(mode="w", encoding="utf-8", dir=out.parent, delete=False) as staging:
            nbformat.write(nb, staging)
        Path(staging.name).replace(out)
    with client.setup_kernel():
        for i, cell in enumerate(nb.cells):
            print(f"{path.name} {i+1}/{len(nb.cells)} {cell.source.splitlines()[0][:78]}", flush=True)
            try:
                client.execute_cell(cell, i)
            finally:
                checkpoint()

if __name__ == "__main__":
    main()
