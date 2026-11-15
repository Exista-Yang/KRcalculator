# INTEG verification report (branch `integrate` @ fe824d7)

Date: 2026-11-15. Machine: 4-core Xeon 2.1 GHz. Result: **no fixes needed**; the merge compiles and passes.

## 1. Build and tests
- `cargo build --release`: ok (warnings only, none new from the merge).
- `cargo build --release --features dev`: ok.
- `KRHOM_THREADS=4 cargo test --release`: 32 test binaries, **320 passed, 0 failed, 43 ignored**.

## 2. Libraries and smoke test
- `tools/fetch_library.sh --root <fresh> l222_n4 l2222_n4`: partial fetch of data commit ee548d2
  (40456 KiB of git objects, 63 s); all 4 library/prime pairs installed and verified.
- `tools/smoke.sh --lib-dir <root> --require-libs`: all 4 checks ok (6 11 / 1 7 22 / 1 3 9 / 1 5 25), 23 s.

## 3. Feature interplay
- (a) `--mark` + `--checkpoint`: braid `1 1 2 -1 2` (3 strands, 2 components), colors 2, `--reduced`, N = 2,3,4.
  `--mark 1` conjugates to `1 2 -1 2 1`; the manifest stores the conjugated braid (`braid=1 2 -1 2 1`),
  so it differs from the `--mark 0` manifest (`braid=1 1 2 -1 2`). Reusing the `--mark 1` directory with
  `--mark 0` is refused (exit 1, "braid: directory has ..."). The checkpointed output equals the plain
  `--reduced --mark 1` run (ranks 1 5 18). The mark is captured through the braid, so no manifest change was made.
- (b) `--knot 4_1 --colors 2 --reduced --checkpoint D2`, killed with `timeout` at 7 s, 9 s and 14 s, then rerun:
  every resumed output equals the plain run (ranks 1 5 25). At 14 s the f31 prime was complete and was
  resumed (36/36 points already done). Kills inside a prime's scan restart that prime's scan, because
  the 36 points of this word are written as one checkpoint file.
- (c) `--knot 6_1 --colors 2 --reduced --estimate`: exit 0 after 491 s. Prints 216 box points per prime,
  about 1h08m wall for 2 primes at 4 threads, peak about 9.6 GB. Note: the estimate itself takes about
  8 minutes here (library load about 3 min per prime). The README does not say this.
- (d) `--file` with 2 braids (aaa, aBaB), N = 2,3, `--latex`: exit 0, 4 tabulars, one section per braid.
- (e) `krhom 33 --braid aaa` and `krhom 36 --braid aaa --colors 3`: both exit 1 with the limit message.

## 4. README
- All 12 relative links/anchors and every `docs/ tools/ data/ examples/ ...` path exist.
- Ran: `2,3,4 --knot 8_19` (8 19 30), `2,3,4 --braid aaa --colors 2` (1 7 22), `2,3 --braid aBaB` (6 11),
  `2,3,4 --braid "1 1 1" --strands 2 --colors 2 --reduced` (1 3 9), `2,3 --link 'L2a1{1}'` (4 9).
  `2,3,4 --knot 4_1 --colors 2 --estimate`: exit 0 (about 57 s for 2 primes). 6_1: see 3(c).
  `--checkpoint` parsing: see 3(a)/(b).

## 5. CI
- `.github/workflows/ci.yml` parses as YAML (job `build-test`). Its build, test (with skips), dev build and
  `smoke.sh --no-libs` steps match commands that pass locally.
