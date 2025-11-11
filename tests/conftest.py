import sys
from pathlib import Path

import matplotlib

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def pytest_sessionstart(session):  # noqa: D401
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 120,
            "font.size": 10,
            "axes.grid": True,
            "axes.facecolor": "white",
        }
    )
