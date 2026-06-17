# CHESS-QC — serverless web front-end (Pyodide)

A static page that runs the **same** `chessqc_*` application Python in the browser via
**Pyodide** (no backend). It reads the same contract (`APP_META`/`INPUTS`/`OUTPUTS`/
`compute`) as the desktop driver and renders the same calculator (inputs left,
Compute, results right; SI/US toggle; canvas plot of the profile arrays; data table).

## Files
- `index.html`, **landing page**: hero, highlights, and a functional-areas overview;
  links to the applications and documentation.
- `apps.html`, **application selection**: lists applications by area with a search
  filter and a fidelity legend; each card links to `calc.html?app=<id>`.
- `calc.html`, **the calculator**; reads `?app=<id>`, with "← Applications" and "Docs" links.
- `docs.html` + `docs.js`, **in-site documentation viewer**: renders the Markdown in
  `../../docs/*.md` with marked (loaded from a CDN), themed by `style.css`.
- `apps.js`, shared **application manifest** (`CHESSQC_APPS` / `CHESSQC_AREAS`) used by
  the landing, selection page, and the driver; add one entry per implemented application.
- `theme.js`, applies mode (light/dark), vibe, badge, and **palette** (Original/Vibrant,
  launcher-set) from `_prefs.js`; exposes `CHESSQCUI` (incl. `wireToggle`).
- `driver.js`, boots Pyodide, loads numpy, fetches the app `.py` + `bridge.py`,
  reads the contract, builds the form, runs `compute()`, renders results + canvas plot.
- `bridge.py`, pure-Python `contract(mod)` / `run(mod, json)` (serialize contract,
  run compute → JSON). Runs identically under Pyodide and CPython (so it is
  unit-testable without a browser, see Verification).
- `style.css`, the shared stylesheet (WaveMaker tokens; mirrors the desktop look).

## Run
Easiest, the top-level launcher (starts the server + opens the hub for you):
```
python run_chess_quick_compute.py --web        # or just run it and pick "Browser"
python run_chess_quick_compute.py --web --app 2-1   # jump straight to one app
```
Or manually (Pyodide needs `fetch()`, so serve over http, not `file://`). From the
**repo root**:
```
python -m http.server 8000
# landing:    http://localhost:8000/chessqc_web/index.html
# apps:       http://localhost:8000/chessqc_web/apps.html
# direct app: http://localhost:8000/chessqc_web/calc.html?app=2-1
```
`calc.html` fetches the application from `../applications/...` (relative to the page),
so the server must be rooted at the repo root.

First load pulls Pyodide + numpy (~10–30 MB) from the jsDelivr CDN, then it's cached.

## Adding applications
Extend the `APPS` manifest in `driver.js`:
```js
const APPS = { "2-1": "../applications/chessqc_2_1_linear_wave_theory.py", /* "2-4": "...", */ };
```
Any application exposing the contract works with no other changes. Apps needing
`scipy` (e.g. cnoidal/Fenton): add `await py.loadPackage("scipy")` in `boot()`.

## Fully offline (no CDN)
1. Download a Pyodide release (`pyodide-<ver>.tar.bz2`) + the `numpy` wheel into `./vendor/`.
2. Point the loader at it: `<script src="vendor/pyodide.js">` and
   `loadPyodide({ indexURL: "vendor/" })`.
3. Bundle the app `.py` + `bridge.py` (already same-origin). The page then runs with
   no network, open it from a share/USB via any static server.

## Verification
- **Python side (done, browser-free):** `bridge.py` was validated under CPython
  against the 2-1 contract, `contract()` emits the right meta/inputs/outputs JSON,
  `run()` returns L=70.88 m, η=1 m, U_r=10.05, 201-point profiles, and the error path
  returns the validation message. This is the exact code Pyodide executes.
- **Browser side (run locally / CI):** load the page per "Run" above; confirm the
  calculator builds, Compute fills the values + plot + table, and the SI/US toggle
  converts (L 70.883 m ↔ 232.56 ft). For CI, a Playwright script can load the page,
  wait for boot, click Compute, and assert the rendered values match the headless
  `compute()` (the same check the desktop driver's smoke test performs).
