#!/usr/bin/env bash
# smoke.sh: quick check that a krhom build gives known values (a few seconds; ~25 s more with the
# library case). Exit 0 if every check passes, 1 on usage errors, 2 on a mismatch or failure.
#
# usage: tools/smoke.sh [--bin PATH] [--lib-dir DIR] [--no-libs | --require-libs]
#   --bin PATH      krhom binary (default $KRHOM_BIN, else target/release/krhom; if that is missing,
#                   `cargo build --release --bin krhom` builds only the calculator)
#   --lib-dir DIR   library root (default $KRHOM_LIB_DIR, else ${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib)
#   --no-libs       skip the library case even if l222_n4 is installed
#   --require-libs  fail (instead of skip) if l222_n4 is not installed
#
# Checks (total rank per N, in order):
#   krhom 2,3 --braid aBaB                                       6 11
#   krhom 2,3,4 --braid aaa --colors 2                           1 7 22
#   krhom 2,3,4 --braid aaa --colors 2 --reduced                 1 3 9
#   krhom 2,3,4 --braid "1 -2 1 -2" --colors 2 --reduced         1 5 25   (library l222_n4 only)
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$here/.." && pwd)"
bin="${KRHOM_BIN:-}"
lib="${KRHOM_LIB_DIR:-${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib}"
libs=auto
usage() { sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 1; }
while [ $# -gt 0 ]; do
  case "$1" in
    --bin) bin="$2"; shift 2 ;;
    --lib-dir) lib="$2"; shift 2 ;;
    --no-libs) libs=no; shift ;;
    --require-libs) libs=yes; shift ;;
    -h|--help) usage ;;
    *) echo "smoke: unknown argument $1" >&2; usage ;;
  esac
done
if [ -z "$bin" ]; then
  bin="$repo_root/target/release/krhom"
  if [ ! -x "$bin" ]; then
    echo "smoke: building the calculator (cargo build --release --bin krhom) ..." >&2
    (cd "$repo_root" && cargo build --release --bin krhom) >&2 || { echo "smoke: build failed" >&2; exit 2; }
  fi
fi
[ -x "$bin" ] || { echo "smoke: no executable krhom at $bin" >&2; exit 2; }

fail=0
t0=$(date +%s)
# check NAME "EXPECTED RANKS" ARGS...: run krhom, compare the "total rank" lines in order
check() {
  local name="$1" want="$2"; shift 2
  local out rc=0 got s=$(date +%s)
  out="$("$bin" "$@" 2>/dev/null)" || rc=$?
  got="$(awk '/^  total rank: /{printf "%s%s", sep, $3; sep=" "}' <<<"$out")"
  if [ "$rc" = 0 ] && [ "$got" = "$want" ]; then
    echo "ok    $name: ranks $got ($(( $(date +%s) - s )) s)"
  else
    echo "FAIL  $name: exit $rc, ranks '$got', expected '$want'  (krhom $*)"
    fail=1
  fi
}
check "4_1 uncolored, N=2,3" "6 11" 2,3 --braid aBaB
check "trefoil color 2, N=2,3,4" "1 7 22" 2,3,4 --braid aaa --colors 2
check "trefoil color 2 reduced, N=2,3,4" "1 3 9" 2,3,4 --braid aaa --colors 2 --reduced

if [ "$libs" = no ]; then
  echo "skip  4_1 color 2 reduced on 3 strands (library l222_n4): --no-libs"
elif "$repo_root/tools/fetch_library.sh" --root "$lib" --verify-only l222_n4 >/dev/null 2>&1; then
  check "4_1 color 2 reduced, 3 strands (library l222_n4)" "1 5 25" 2,3,4 --braid "1 -2 1 -2" --colors 2 --reduced --lib-dir "$lib"
elif [ "$libs" = yes ]; then
  echo "FAIL  library l222_n4 is not installed (or does not verify) under $lib: run tools/fetch_library.sh l222_n4"
  fail=1
else
  echo "skip  4_1 color 2 reduced on 3 strands: library l222_n4 not installed under $lib (tools/fetch_library.sh l222_n4 installs it)"
fi

if [ "$fail" = 0 ]; then echo "smoke: all checks passed ($(( $(date +%s) - t0 )) s)"; exit 0; fi
echo "smoke: FAILED" >&2
exit 2
