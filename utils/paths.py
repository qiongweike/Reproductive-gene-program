"""Project directories. Run scripts from the repository root; set SPG_ROOT to keep data and results elsewhere."""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("SPG_ROOT", REPO))
DATA = ROOT / "data"
RAW = DATA / "raw"
METADATA = DATA / "metadata"
PROCESSED = DATA / "processed"
EXTERNAL = DATA / "external"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
VENDOR = ROOT / "vendor"
RESOURCES = REPO / "resources"


def out_dir(*parts):
    """Return RESULTS/parts, creating it if needed."""
    d = RESULTS.joinpath(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d
