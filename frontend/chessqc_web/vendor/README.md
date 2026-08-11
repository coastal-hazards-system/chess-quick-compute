# Vendored front-end libraries

Third-party code kept in the repository rather than fetched from a CDN, so the
calculator keeps working offline, under a strict content-security policy, and
from any static host.

## hyparquet

- `hyparquet.esm.js` — hyparquet 1.28.1, bundled ESM build (jsDelivr/Rollup).
- `hyparquet.LICENSE` — MIT, Hyperparam.
- Upstream: https://github.com/hyparam/hyparquet

Reads the Parquet station records in `data/` and any Parquet file the user
uploads: the column data and the `pystorm` key of the file footer, which carries
the station id, product, units and **vertical datum** that the calculator
displays and passes down the workflow.

Snappy decompression is built into hyparquet itself, so the companion
`hyparquet-compressors` package is not needed; PyStorm writes Parquet v2.6 with
the SNAPPY codec.

To update: fetch the same path at the new version, keep the license file beside
it, and re-run `tests/test_browser.py`.
