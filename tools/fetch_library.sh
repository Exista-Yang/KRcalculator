#!/usr/bin/env bash
# fetch_library.sh: install (or verify) the prebuilt layer-A libraries the krhom calculator needs
# for colored braids on 3 and 4 strands (uniform color 2, N <= 4).
#
# usage: tools/fetch_library.sh [--root DIR] [--git-dir DIR | --repo URL] [--force] [--verify-only] ID[/field] ...
#   ID      l222_n4 (3 strands) or l2222_n4 (4 strands); /f31 or /f32b selects one prime (default both)
#   --root  install root (default $KRHOM_LIB_DIR, else ${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib)
#   --git-dir DIR  git clone that has (or can fetch) the data commit (default: this checkout)
#   --repo URL     fetch the data commit from URL into a temporary bare repository instead
#   --force        reinstall even if an installed copy verifies (the old copy is moved aside, then removed)
#   --verify-only  only check the installed copy: LIBRARY.txt against the registry, every listed
#                  file present with its sha256 (extra files, e.g. templates written later, are fine)
#
# The registry is data/libraries.tsv (shared with src/calc.rs). Every chunk is checked against the
# registry sha256 and the MANIFEST sha256 of the data commit, its tar listing is checked (regular
# files under layer_a/ and tmpl/ only, counts = MANIFEST), layers are extracted in registry
# (restore) order, and LIBRARY.txt records a per-file sha256 list. Exit 0 on success, 1 on usage
# errors, 2 on a failed download or verification.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$here/.." && pwd)"
registry="$repo_root/data/libraries.tsv"

root="${KRHOM_LIB_DIR:-${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib}"
git_dir="$repo_root"
repo_url=""
force=0
verify_only=0
ids=()
die() { echo "fetch_library: $*" >&2; exit 2; }
usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 1; }
while [ $# -gt 0 ]; do
  case "$1" in
    --root) root="$2"; shift 2 ;;
    --git-dir) git_dir="$2"; shift 2 ;;
    --repo) repo_url="$2"; shift 2 ;;
    --force) force=1; shift ;;
    --verify-only) verify_only=1; shift ;;
    -h|--help) usage ;;
    -*) echo "fetch_library: unknown option $1" >&2; usage ;;
    *) ids+=("$1"); shift ;;
  esac
done
[ ${#ids[@]} -gt 0 ] || usage
[ -f "$registry" ] || die "registry $registry not found"

if command -v sha256sum >/dev/null 2>&1; then sha() { sha256sum "$@"; }
elif command -v shasum >/dev/null 2>&1; then sha() { shasum -a 256 "$@"; }
else die "need sha256sum or shasum"; fi
sha_of() { sha "$1" | awk '{print $1}'; }
field_name() { case "$1" in f31) echo F_2147483647 ;; f32b) echo F_4294967291 ;; *) die "unknown field $1" ;; esac; }
code_commit="$(git -C "$repo_root" rev-parse --short HEAD 2>/dev/null || echo unknown)"

# registry rows matching an id (and optional field)
rows_for() {
  local id="${1%%/*}" fld=""
  [ "$1" != "$id" ] && fld="${1#*/}"
  awk -F'\t' -v id="$id" -v f="$fld" '!/^#/ && NF>=10 && $1==id && (f=="" || $2==f)' "$registry"
}

# --- verification of an installed copy ------------------------------------------------------
# verify DEST ROW: 0 if LIBRARY.txt matches the registry row and every listed file has its hash
verify_installed() {
  local dest="$1" row="$2" quiet="${3:-0}"
  local id fld lam n prov commit costs status stl layers
  IFS=$'\t' read -r id fld lam n prov commit costs status stl layers <<<"$row"
  local lib="$dest/LIBRARY.txt"
  [ -f "$lib" ] || { [ "$quiet" = 1 ] || echo "  $id/$fld: no LIBRARY.txt under $dest" >&2; return 1; }
  local k v
  for kv in "format:krhom-library v1" "id:$id" "lambda:$lam" "n:$n" "field:$(field_name "$fld")" "provider:$prov" "data_commit:$commit" "layers:$layers"; do
    k="${kv%%:*}"; v="${kv#*:}"
    if [ "$(grep -m1 "^$k: " "$lib" | sed "s/^$k: //")" != "$v" ]; then
      [ "$quiet" = 1 ] || echo "  $id/$fld: LIBRARY.txt $k differs from the registry (expected '$v')" >&2
      return 1
    fi
  done
  local cs; cs="$(grep -m1 '^costs: ' "$lib" | awk '{print $3}')"
  [ -f "$dest/costs.txt" ] && [ "$(sha_of "$dest/costs.txt")" = "$cs" ] || { [ "$quiet" = 1 ] || echo "  $id/$fld: costs.txt missing or changed" >&2; return 1; }
  head -1 "$dest/costs.txt" | grep -q '^krhom-costs v2 alg=hs-kl1' || { [ "$quiet" = 1 ] || echo "  $id/$fld: costs.txt has the wrong format" >&2; return 1; }
  local nlisted; nlisted="$(sed -n '/^files:$/,$p' "$lib" | tail -n +2 | wc -l | tr -d ' ')"
  local want; want=$(( $(grep -m1 '^files_layer_a: ' "$lib" | awk '{print $2}') + $(grep -m1 '^files_tmpl: ' "$lib" | awk '{print $2}') ))
  [ "$nlisted" = "$want" ] || { [ "$quiet" = 1 ] || echo "  $id/$fld: LIBRARY.txt lists $nlisted files, header says $want" >&2; return 1; }
  local bad
  bad="$(cd "$dest" && sed -n '/^files:$/,$p' LIBRARY.txt | tail -n +2 | awk '{print $2"  "$3}' | sha -c --quiet 2>&1 || true)"
  if [ -n "$bad" ]; then
    [ "$quiet" = 1 ] || { echo "  $id/$fld: $(echo "$bad" | grep -c . ) listed file(s) missing or changed, e.g.:" >&2; echo "$bad" | head -3 | sed 's/^/    /' >&2; }
    return 1
  fi
  return 0
}

# --- git transport ----------------------------------------------------------------------------
tmp_repo=""
cleanup() { [ -n "$tmp_repo" ] && rm -rf "$tmp_repo"; [ -n "${staging:-}" ] && rm -rf "$staging"; true; }
trap cleanup EXIT
ensure_commit() {
  local commit="$1"
  if [ -n "$repo_url" ] && [ -z "$tmp_repo" ]; then
    tmp_repo="$(mktemp -d)"
    git init -q --bare "$tmp_repo"
    git_dir="$tmp_repo"
  fi
  if git -C "$git_dir" cat-file -e "$commit^{commit}" 2>/dev/null; then return 0; fi
  local remote="${repo_url:-origin}"
  echo "fetching data commit ${commit:0:12} from $remote ..." >&2
  # by SHA (review K11: a branch fetch loses the commit once the branch moves); fall back to the
  # data-gl6 branch with full history if the server refuses SHA fetches
  git -C "$git_dir" fetch -q --no-tags "$remote" "$commit" 2>/dev/null \
    || git -C "$git_dir" fetch -q --no-tags "$remote" "refs/heads/data-gl6:refs/remotes/krhom-fetch/data-gl6" \
    || die "cannot fetch data commit $commit from $remote"
  git -C "$git_dir" cat-file -e "$commit^{commit}" 2>/dev/null || die "data commit $commit not found in $remote (data-gl6 history)"
}

# --- install one registry row -----------------------------------------------------------------
install_row() {
  local row="$1"
  local id fld lam n prov commit costs status stl layers
  IFS=$'\t' read -r id fld lam n prov commit costs status stl layers <<<"$row"
  local dest="$root/$id/$fld"
  if [ "$verify_only" = 1 ]; then
    if verify_installed "$dest" "$row"; then echo "$id/$fld: verified ($dest)"; return 0; else echo "$id/$fld: NOT OK ($dest): run tools/fetch_library.sh --force $id/$fld" >&2; return 2; fi
  fi
  if [ "$force" = 0 ] && [ -d "$dest" ]; then
    if verify_installed "$dest" "$row" 1; then echo "$id/$fld: already installed ($dest)"; return 0; fi
    echo "$id/$fld: an installed copy under $dest does not verify; reinstalling" >&2
  fi
  ensure_commit "$commit"
  mkdir -p "$root"
  staging="$(mktemp -d "$root/.staging-$id-$fld-XXXXXX")"
  mkdir -p "$staging/chunks" "$staging/lib/layer_a" "$staging/lib/tmpl" "$staging/lib/MANIFESTS"
  local i=0 entry path want dir manifest base got na nt man_la man_tm
  IFS=';' read -r -a entries <<<"$layers"
  for entry in "${entries[@]}"; do
    path="${entry%@*}"; want="${entry##*@}"; dir="$(dirname "$path")"
    manifest="$staging/MANIFEST.$i"
    git -C "$git_dir" show "$commit:$dir/MANIFEST.txt" >"$manifest" || die "no MANIFEST for $dir in $commit"
    base="$(basename "$path")"
    git -C "$git_dir" show "$commit:$path" >"$staging/chunks/$i" || die "cannot read $path from $commit"
    got="$(sha_of "$staging/chunks/$i")"
    [ "$got" = "$want" ] || die "$path: sha256 $got, registry says $want"
    grep -q "^$got  $base\$" "$manifest" || die "$path: sha256 $got is not the MANIFEST value"
    # tar listing: regular files and directories under layer_a/ and tmpl/ only
    local listing; listing="$(tar -tvzf "$staging/chunks/$i")" || die "$path: not a gzip tar"
    local names; names="$(tar -tzf "$staging/chunks/$i")"
    if echo "$listing" | awk '{print substr($1,1,1)}' | grep -qv '^[-d]$'; then die "$path: tar has entries that are not regular files or directories"; fi
    if echo "$names" | grep -Evq '^(\./)?((layer_a|tmpl)/?|layer_a/a-[0-9a-f]{16}\.txt|tmpl/tmpl-[0-9a-f]{16}\.txt|\./?)$'; then
      die "$path: unexpected tar entry: $(echo "$names" | grep -Ev '^(\./)?((layer_a|tmpl)/?|layer_a/a-[0-9a-f]{16}\.txt|tmpl/tmpl-[0-9a-f]{16}\.txt|\./?)$' | head -1)"
    fi
    na="$(echo "$names" | grep -Ec 'layer_a/a-[0-9a-f]{16}\.txt$' || true)"
    nt="$(echo "$names" | grep -Ec 'tmpl/tmpl-[0-9a-f]{16}\.txt$' || true)"
    # MANIFEST "files:" line in either word order ("layer_a 233, tmpl 471" / "249 layer_a, 233 tmpl")
    local fl; fl="$(grep -m1 '^files' "$manifest")"
    man_la="$(echo "$fl" | grep -oE '(layer_a [0-9]+|[0-9]+ layer_a)' | grep -oE '[0-9]+')"
    man_tm="$(echo "$fl" | grep -oE '(tmpl [0-9]+|[0-9]+ tmpl)' | grep -oE '[0-9]+')"
    [ "$na" = "$man_la" ] && [ "$nt" = "$man_tm" ] || die "$path: tar has $na layer_a / $nt tmpl files, MANIFEST says $man_la / $man_tm"
    echo "$names" | grep -E '(layer_a/a-|tmpl/tmpl-)' | sed 's#^\./##' | sed "s#^#$i #" >>"$staging/owners"
    tar -C "$staging/lib" -xzf "$staging/chunks/$i" || die "$path: extraction failed"
    cp "$manifest" "$staging/lib/MANIFESTS/$(echo "$dir" | tr '/' '_').txt"
    echo "  $id/$fld layer $i: $path ok (sha256 ${got:0:12}, $na layer_a + $nt tmpl)" >&2
    i=$((i + 1))
  done
  git -C "$git_dir" show "$commit:$costs" >"$staging/lib/costs.txt" || die "cannot read $costs"
  head -1 "$staging/lib/costs.txt" | grep -q '^krhom-costs v2 alg=hs-kl1' || die "$costs: wrong format"
  local status_line="complete"
  if [ "$status" != "complete" ]; then
    [ -f "$repo_root/$status" ] || die "status file $status missing from this checkout"
    cp "$repo_root/$status" "$staging/lib/status.txt"
    status_line="status.txt"
  fi
  # last writer of every file (registry order: later layers overwrite earlier names)
  awk '{own[$2]=$1} END {for (f in own) print own[f], f}' "$staging/owners" | sort -k2 >"$staging/final"
  local fa ft
  fa="$(grep -c ' layer_a/' "$staging/final" || true)"; ft="$(grep -c ' tmpl/' "$staging/final" || true)"
  [ "$fa" = "$(find "$staging/lib/layer_a" -type f | wc -l | tr -d ' ')" ] || die "layer_a count mismatch after extraction"
  [ "$ft" = "$(find "$staging/lib/tmpl" -type f | wc -l | tr -d ' ')" ] || die "tmpl count mismatch after extraction"
  {
    echo "format: krhom-library v1"
    echo "id: $id"
    echo "lambda: $lam"
    echo "n: $n"
    echo "field: $(field_name "$fld")"
    echo "provider: $prov"
    echo "data_commit: $commit"
    echo "layers: $layers"
    echo "stamp_layers: $stl"
    echo "costs: costs.txt $(sha_of "$staging/lib/costs.txt")"
    echo "status: $status_line"
    echo "files_layer_a: $fa"
    echo "files_tmpl: $ft"
    echo "installed: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "fetched_by: tools/fetch_library.sh @ $code_commit"
    echo "files:"
    (cd "$staging/lib" && while read -r own f; do echo "$own $(sha_of "$f") $f"; done <"$staging/final")
  } >"$staging/lib/LIBRARY.txt"
  local aside=""
  if [ -d "$dest" ]; then aside="$dest.old-$$"; mv "$dest" "$aside"; fi
  mkdir -p "$(dirname "$dest")"
  mv "$staging/lib" "$dest"
  [ -n "$aside" ] && rm -rf "$aside"
  rm -rf "$staging"; staging=""
  verify_installed "$dest" "$row" || die "$id/$fld: installed copy does not verify"
  echo "$id/$fld: installed ($dest: $fa layer_a, $ft tmpl files, data ${commit:0:12})"
}

rc=0
for want in "${ids[@]}"; do
  rows="$(rows_for "$want")"
  [ -n "$rows" ] || { echo "fetch_library: unknown library $want (known: $(awk -F'\t' '!/^#/ && NF>=10 {print $1"/"$2}' "$registry" | tr '\n' ' '))" >&2; exit 1; }
  while IFS= read -r row; do
    install_row "$row" || rc=$?
  done <<<"$rows"
done
exit $rc
