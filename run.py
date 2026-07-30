"""Single entry point for the charge-transfer cluster calculations.

Set the two variables below and press F5 in Spyder, or pass them on the
command line:

    python run.py                                  # uses the settings below
    python run.py inputs/Fe_Ba2FeSi2O7.py          # that material, same mode
    python run.py inputs/Fe_Ba2FeSi2O7.py xps      # that material, XPS only
    python run.py inputs/NiO.py both               # XAS and XPS in one go

One input file describes one material. Parameters the two spectroscopies share
-- the cluster itself -- are written once; anything that differs carries an
`xas_` or `xps_` prefix. Running `both` therefore guarantees the two spectra
come from exactly the same cluster, which is the point of fitting them
together.

`run_xas.py` and `run_xps.py` still work on their own and read the same input
files; this script only saves you from editing `INPUT_FILE` in two places.
"""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# What to run.

INPUT_FILE = ROOT / "inputs" / "Fe_Ba2FeSi2O7.py"
SPECTROSCOPY = "xas"  # "xas", "xps", or "both"


_MODES = ("xas", "xps", "both")


def main(argv: list[str] | None = None) -> None:
    input_file, spectroscopy = _resolve_arguments(argv)
    modes = ("xas", "xps") if spectroscopy == "both" else (spectroscopy,)

    for index, mode in enumerate(modes):
        if index:
            print()
        print("=" * 72)
        print(f"  {mode.upper()}   {Path(input_file).name}")
        print("=" * 72)
        _run_one(mode, input_file)


def _run_one(mode: str, input_file) -> None:
    if mode == "xas":
        import run_xas as runner
    else:
        import run_xps as runner
    runner.main(input_file)


def _resolve_arguments(argv: list[str] | None) -> tuple[object, str]:
    arguments = list(sys.argv[1:] if argv is None else argv)
    input_file = INPUT_FILE
    spectroscopy = SPECTROSCOPY

    for argument in arguments:
        lowered = argument.strip().lower()
        if lowered in {"-h", "--help", "help"}:
            print(_usage())
            raise SystemExit(0)
        if lowered in _MODES:
            spectroscopy = lowered
        else:
            path = Path(argument)
            if not path.is_absolute():
                path = ROOT / path
            if not path.exists():
                raise SystemExit(f"input file not found: {path}")
            input_file = path

    if spectroscopy not in _MODES:
        raise SystemExit(f"spectroscopy must be one of {_MODES}, got {spectroscopy!r}")
    return input_file, spectroscopy


def _usage() -> str:
    available = sorted(path.name for path in (ROOT / "inputs").glob("*.py"))
    return "\n".join(
        [
            "usage: python run.py [input file] [xas | xps | both]",
            "",
            "Arguments may be given in either order and both are optional;",
            "whatever is omitted falls back to the settings at the top of this file",
            f"(currently {Path(INPUT_FILE).name}, {SPECTROSCOPY}).",
            "",
            "available input files: " + (", ".join(available) or "none"),
        ]
    )


if __name__ == "__main__":
    main()
