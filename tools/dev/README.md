# Developer tools

Scripts for library builds, sharded box runs and their analysis. They need the developer
binaries (`cargo build --release --features dev`). The user tools are one level up
(`tools/fetch_library.sh`, `randcheck.py`, `batch_check.py`, `wang.py`). Scripts that only
served the archived runs are on branch `archive-pre-restructure` (`tools/`, `tools/c2/`).

| script | purpose |
|---|---|
| `box.sh`, `finish.sh` | sharded 4-strand (2,2,2,2) box run with `sbox box` / `sbox finish` (how 6_1, 7_6, 7_7 in `results/colored_gl4/` were computed); paths from `CACHE`, `COSTS`, `OUTD` |
| `persist_layer_a.sh`, `persist_c2.sh` | pack a layer-A / tmpl cache into ≤ 45 MB tar.gz chunks with a MANIFEST for the `data-gl6` branch |
| `libbuild_report.py` | summarize `libbuild` run files (JSON lines) |
| `trace_windows.py` | refused Search windows from `KRHOM_TRACE` logs → `libbuild calib --windows` lines |
| `coverx.py` | X distribution of the library jobs each knot needs (`libbuild cover --list` + plan JSON) and an expected-failure projection |
| `xdist.py`, `xfail.py` | X-flow distribution of `libbuild` plans; outcome and wall per X bucket of `libbuild` runs |
| `k4m_status.py` | job status lines (`JOB ok\|fail`) from `libbuild` / runjobs JSON lines (the format of `data/libstatus/`) |
| `same.py` | are all closures of two `sbox` / `xscan` JSON records equal? |
| `xscan_vs_wang.py` | `xscan` 2-strand unreduced closures against `tests/data/wang_ref.txt` |
| `ntrunc_t1/` | separate scratch crate of the N-truncation study (`cargo run --release --manifest-path tools/dev/ntrunc_t1/Cargo.toml`); `data/` holds the 3- and 4-strand knot braid lists |
