# Development documents

Design notes, specifications, run reports and independent reviews written while the calculator
was built. The user documentation is [README.md](../../README.md) and
[docs/CALCULATOR.md](../CALCULATOR.md); recorded results are in [results/](../../results/README.md).

Paths of the form `runs_*/...` in these documents (scripts, queues, logs, intermediate JSON)
refer to the development directories that are kept on branch `archive-pre-restructure`
(`git show archive-pre-restructure:runs_c2/...`). Source paths such as `src/sindec.rs` refer to
the flat layout before the restructure; the module now lives in a group directory (for example
`src/library/sindec.rs`); the old → new table is in [RESTRUCTURE_report.md](RESTRUCTURE_report.md).

| path | content |
|---|---|
| `PLAN.md` | goals, model, fixed conventions (gradings, crossing signs), milestone history |
| `RESTRUCTURE_report.md` | the repository restructure: path table, deleted files, stamp check, validation |
| `colored_design/SYNTHESIS.md` | start here: the chosen colored gl_N design and plan |
| `colored_design/SPEC_M4_M6a.md`, `M4a_design.md`, `M4b_design.md`, `M6_0_design.md`, `HANDOFF_M4b.md` | specifications and design notes of the Hom provider (M4), the exact colored engine (M6.0) and the localized engine (M6a) |
| `colored_design/NTRUNC_STUDY.md` | N-truncation study (tool: `tools/dev/ntrunc_t1`) |
| `colored_design/gl6/` | (3,3,3) gl_6 work: specification rev1, milestone reports M5a–M6a-2a, MBOX, the 2-strand σ^9 data (`m2s/`) |
| `colored_design/m4c/` | light-leaves backend (M4c): plan, reports A1–B4w, data |
| `colored_design/b5/` | B5 census and analysis |
| `colored_design/c2/` | color-2 gl_4 sessions C2-1 … C2-7 (3- and 4-strand libraries and results) |
| `colored_design/calc/` | the calculator (CALC-1 … CALC-3): specification, decisions, reviews, reports |
| `colored_design/review_*` | independent review waves |
