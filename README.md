# CHESS-QC
*Coastal Hazards, Engineering, and Structures System (CHESS) — Quick Compute (QC)*

A self-contained toolkit of **35+ coastal-engineering calculators**: an open-source
reimplementation of the U.S. Army Corps of Engineers **ACES** (Automated Coastal
Engineering System) computational library, runnable in the browser with no
install, as a native desktop application, or as standalone Python scripts.

Each application reproduces an ACES routine across wave prediction and theory, wave
transformation, structural design, runup/transmission/overtopping, littoral processes,
inlet and harbor design, and storm surge. Every tool is validated against the ACES
*User's Guide* worked examples or, where none exist, against analytic limits and the
primary literature.

---

## Highlights

- **Serverless web calculator:** runs entirely in the browser via [Pyodide](https://pyodide.org), with no backend and no install. The web driver executes the same Python application files the desktop and command-line front-ends use.
- **Native desktop application:** a PySide6 hub with per-application calculators.
- **Standalone scripts:** every application is a single, dependency-light `.py` (numpy plus the standard library) that runs its own validation and prints a worked example.
- **SI and US units** throughout, converted at the edge (SI internally).
- **Fidelity labels:** each application is classified *exact*, *standard*, or *provisional*, with the rationale stated in its header and in the validation report.
- **Documentation that stays in sync:** the User Manual and Validation Report are generated directly from the application contracts.

## Functional areas

| Area | Applications | Count |
|---|---|---|
| Wave Prediction | 1-1 to 1-6 | 6 |
| Wave Theory | 2-1 to 2-5 | 5 |
| Wave Transformation | 3-1 to 3-4 | 4 |
| Structural Design | 4-1 to 4-4 | 4 |
| Wave Runup, Transmission and Overtopping | 5-1 to 5-5 | 5 |
| Littoral Processes | 6-1 to 6-5 | 5 |
| Inlet Processes | 7-1, 7-2 | 2 |
| Harbor Design | 8-1, 8-2, 8-3 | 3 |
| Storm Surge | 9-1 | 1 |
| Miscellaneous Routines | M-1 | 1 |

See [`docs/USER_MANUAL.md`](docs/USER_MANUAL.md) for the full list with inputs and outputs.

## Quick start

Requires **Python 3.10+**.

```bash
# core dependency (used by the applications)
pip install numpy
```

**Launch the web calculator** (opens a hub page in your browser):

```bash
python run_chess_quick_compute.py --web
# jump straight to one application:
python run_chess_quick_compute.py --web --app 2-1
```

> The browser version fetches Pyodide and numpy from a CDN on first load (this needs
> internet once; afterwards it is cached).

**Launch the native desktop application** (requires PySide6):

```bash
pip install PySide6 numpy
python run_chess_quick_compute.py --desktop
```

Running `python run_chess_quick_compute.py` with no flags prompts you to choose a
front-end.

**Run a single application from the command line.** It executes its self-tests and prints
the worked example:

```bash
python backend/applications/chessqc_2_1_linear_wave_theory.py
```

## Fidelity classification

Every application carries one of three labels (shorthand **I**, **II**, **III** in the docs):

- **I (exact):** every coefficient and relationship is known from the source and the results are validated, with nothing inferred. Empirical formulas qualify when their coefficients are known and verified.
- **II (standard):** a named published method that involves a self-made convention or only partial validation.
- **III (provisional):** a coefficient or relationship had to be inferred, is not recoverable from the public sources, or has no numeric reference example.

Current distribution: **28 exact, 3 standard, 5 provisional.** Each application's full
rationale is in its header and in [`docs/VALIDATION_REPORT.md`](docs/VALIDATION_REPORT.md).

## Validation

The primary reference is the worked "Example Problem" in the ACES *User's Guide*,
reproduced to the stated tolerance. Where ACES gives no numeric example, validation is
analytic (closed-form limits and cross-checks) or against the primary literature.
Documented residuals and caveats are stated per application. Each application also runs its
own self-tests when executed directly, and the per-application reference, method, and
tolerance are tabulated in [`docs/VALIDATION_REPORT.md`](docs/VALIDATION_REPORT.md).

## Project layout

```
backend/applications/   35+ self-contained application modules (chessqc_*.py)
common/                 Python/JS bridge, unit conversion, documentation generators
frontend/chessqc_web/   serverless web calculator (Pyodide)
frontend/chessqc_gui/   PySide6 desktop application
docs/                   User Manual, Technical Reference, Equation Reference,
                        Validation Report, Outdated Applications
run_chess_quick_compute.py   launcher (web or desktop)
```

Each application is **self-contained** (zero sibling imports) and exposes a uniform
contract (`APP_META`, `INPUTS`, `OUTPUTS`, and a `compute()` function) that the web and
desktop front-ends render generically.

## Documentation

- [User Manual](docs/USER_MANUAL.md): per-application inputs, outputs, ranges, and defaults.
- [Technical Reference](docs/TECHNICAL_REFERENCE.md): the physics, assumptions, and caveats behind each tool.
- [Equation Reference](docs/EQUATIONS.md): the equations, with their sources.
- [Validation Report](docs/VALIDATION_REPORT.md): the reference, method, and tolerance for every application.
- [Outdated / Superseded Applications](docs/OUTDATED_APPS.md): tools kept for screening or legacy consistency.

## Disclaimer

CHESS-QC reimplements the methods of ACES. Each application cites ACES and its underlying
published sources (the Shore Protection Manual 1984, the Coastal Engineering Manual, and the
original authors) as the reference for every method. It is intended for engineering
screening, teaching, and research; results should be reviewed by a qualified coastal
engineer before use in design.