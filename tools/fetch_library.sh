#!/usr/bin/env bash
# fetch_library.sh: install (or verify) the prebuilt layer-A libraries the krhom calculator needs
# for colored braids on 3 and 4 strands (uniform color 2, N <= 4).
#
# usage: tools/fetch_library.sh [--root DIR] [--git-dir DIR | --repo URL | --url BASE] [--full-fetch]
#                               [--force] [--verify-only] ID[/field] ...
#   ID      l222_n4 (3 strands) or l2222_n4 (4 strands); /f31 or /f32b selects one prime (default both)
#   --root  install root (default $KRHOM_LIB_DIR, else ${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib)
#   --git-dir DIR  git clone that has (or can fetch from its origin) the data commit (default: this checkout)
#   --repo URL     fetch the data commit from URL (default: the origin of --git-dir)
#   --url BASE     download release assets BASE/<asset> (https://, http:// or file://; curl or wget)
#                  listed in data/release_assets.tsv instead of using git
#   --full-fetch   git: fetch the full history of the data commit (the pre-0.2 method) instead of a
#                  partial fetch (only the needed blobs)
#   --force        reinstall even if an installed copy verifies (the old copy is moved aside, then removed)
#   --verify-only  only check the installed copy: LIBRARY.txt against the registry, every listed
#                  file present with its sha256 (extra files, e.g. templates written later, are fine)
#
set -euo pipefail

# git on a directory that may be a bare repository: with safe.bareRepository=explicit,
# `git -C BARE` refuses implicit discovery, so bare repositories get an explicit --git-dir.
gitd() { local d="$1"; shift
  if [ -f "$d/HEAD" ] && [ -d "$d/objects" ] && [ ! -e "$d/.git" ]; then git --git-dir="$d" "$@"; else git -C "$d" "$@"; fi; }

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$here/.." && pwd)"
registry="$repo_root/data/libraries.tsv"

root="${KRHOM_LIB_DIR:-${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib}"
git_dir="$repo_root"
repo_url=""
url_base=""
full_fetch=0
force=0
verify_only=0
ids=()
die() { echo "fetch_library: $*" >&2; exit 2; }
usage() { sed -n '2,23p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 1; }
while [ $# -gt 0 ]; do
  case "$1" in
    --root) root="$2"; shift 2 ;;
    --git-dir) git_dir="$2"; shift 2 ;;
    --repo) repo_url="$2"; shift 2 ;;
    --url) url_base="${2%/}"; shift 2 ;;
    --full-fetch) full_fetch=1; shift ;;
    --force) force=1; shift ;;
    --verify-only) verify_only=1; shift ;;
    -h|--help) usage ;;
    -*) echo "fetch_library: unknown option $1" >&2; usage ;;
    *) ids+=("$1"); shift ;;
  esac
done
[ ${#ids[@]} -gt 0 ] || usage
if [ -n "$url_base" ] && { [ -n "$repo_url" ] || [ "$full_fetch" = 1 ]; }; then echo "fetch_library: --url excludes --repo and --full-fetch" >&2; exit 1; fi
case "$url_base" in ""|https://*|http://*|file://*) ;; *) echo "fetch_library: --url needs an https://, http:// or file:// base" >&2; exit 1 ;; esac
git_dir0="$git_dir"
assets_tsv="$repo_root/data/release_assets.tsv"
[ -f "$registry" ] || die "registry $registry not found"

if command -v sha256sum >/dev/null 2>&1; then sha_cmd=(sha256sum)
elif command -v shasum >/dev/null 2>&1; then sha_cmd=(shasum -a 256)
else die "need sha256sum or shasum"; fi
sha() { "${sha_cmd[@]}" "$@"; }
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
  # sources of costs.txt and status.txt (review tools #4): the costs path in the data commit must be
  # the registry's (same data_commit + same path = same bytes), and status.txt must be the
  # registry's status file with the checkout's current sha256 (a corrected data/libstatus file
  # makes the install stale)
  if [ "$(grep -m1 '^costs_source: ' "$lib" | sed 's/^costs_source: //')" != "$costs" ]; then
    [ "$quiet" = 1 ] || echo "  $id/$fld: LIBRARY.txt costs_source differs from the registry (expected '$costs')" >&2
    return 1
  fi
  local want_status="complete"
  if [ "$status" != "complete" ]; then
    [ -f "$repo_root/$status" ] || { [ "$quiet" = 1 ] || echo "  $id/$fld: status file $status missing from this checkout" >&2; return 1; }
    want_status="$status $(sha_of "$repo_root/$status")"
    [ -f "$dest/status.txt" ] && [ "$(sha_of "$dest/status.txt")" = "${want_status##* }" ] || { [ "$quiet" = 1 ] || echo "  $id/$fld: status.txt missing or differs from $status in this checkout" >&2; return 1; }
  fi
  if [ "$(grep -m1 '^status_source: ' "$lib" | sed 's/^status_source: //')" != "$want_status" ]; then
    [ "$quiet" = 1 ] || echo "  $id/$fld: LIBRARY.txt status_source differs from the registry/checkout (expected '$want_status')" >&2
    return 1
  fi
  local cs; cs="$(grep -m1 '^costs: ' "$lib" | awk '{print $3}')"
  [ -f "$dest/costs.txt" ] && [ "$(sha_of "$dest/costs.txt")" = "$cs" ] || { [ "$quiet" = 1 ] || echo "  $id/$fld: costs.txt missing or changed" >&2; return 1; }
  head -1 "$dest/costs.txt" | grep -q '^krhom-costs v2 alg=hs-kl1' || { [ "$quiet" = 1 ] || echo "  $id/$fld: costs.txt has the wrong format" >&2; return 1; }
  local nlisted; nlisted="$(sed -n '/^files:$/,$p' "$lib" | tail -n +2 | wc -l | tr -d ' ')"
  local want; want=$(( $(grep -m1 '^files_layer_a: ' "$lib" | awk '{print $2}') + $(grep -m1 '^files_tmpl: ' "$lib" | awk '{print $2}') ))
  [ "$nlisted" = "$want" ] || { [ "$quiet" = 1 ] || echo "  $id/$fld: LIBRARY.txt lists $nlisted files, header says $want" >&2; return 1; }
  local bad
  # Portable check (GNU sha256sum, BSD sha256sum on macOS, perl shasum): hash the listed files
  # and compare with the list; no -c/--quiet flags, whose spelling differs between tools.
  bad="$(cd "$dest" && sed -n '/^files:$/,$p' LIBRARY.txt | tail -n +2 | awk '{print $2, $3}' | {
      want=$(mktemp) && cat > "$want"
      awk '{print $2}' "$want" | while IFS= read -r f; do [ -f "$f" ] || echo "$f: FAILED (missing)"; done
      awk '{print $2}' "$want" | while IFS= read -r f; do [ -f "$f" ] && printf '%s\0' "$f"; done \
        | xargs -0 "${sha_cmd[@]}" 2>/dev/null \
        | awk 'NR==FNR{h[$2]=$1; next} ($2 in h) && h[$2]!=$1 {print $2": FAILED"}' "$want" -
      rm -f "$want"; } || true)"
  if [ -n "$bad" ]; then
    [ "$quiet" = 1 ] || { echo "  $id/$fld: $(echo "$bad" | grep -c "FAILED") listed file(s) missing or changed, e.g.:" >&2; echo "$bad" | head -3 | sed 's/^/    /' >&2; }
    return 1
  fi
  return 0
}

# --- transport: git (partial or full fetch) or release assets ------------------------------
# The bytes come from one of three sources, all checked identically afterwards (chunk sha256 =
# registry = MANIFEST, tar listing, counts):
#   * git, commit already present in --git-dir: read there, nothing is downloaded;
#   * git, partial fetch (default): a temporary bare repository fetches the data commit with
#     --filter=blob:none --depth 1 (trees only, a few KB) and `git show` then fetches exactly the
#     blobs the registry row needs (chunks, MANIFESTs, costs), lazily from the promisor remote;
#     the user's checkout is not touched. If the server refuses partial or SHA fetches, the old
#     method below is used.
#   * git, full fetch (--full-fetch or fallback): the data commit with its whole history, by SHA
#     (review K11: a branch fetch loses the commit once the branch moves), into the checkout (or,
#     with --repo, into a temporary bare repository); last resort: the data-gl6 branch.
#   * --url BASE: one release asset per (library, prime), its sha256 checked against
#     data/release_assets.tsv before extraction, extracted (regular files at the expected data
#     commit paths only) into the staging directory.
tmp_repo=""
fetch_mode=""
cleanup() { [ -n "$tmp_repo" ] && rm -rf "$tmp_repo"; [ -n "${staging:-}" ] && rm -rf "$staging"; true; }
trap cleanup EXIT
new_tmp_repo() {
  [ -n "$tmp_repo" ] && rm -rf "$tmp_repo"
  tmp_repo="$(mktemp -d)" || die "mktemp failed"
  [ -n "$tmp_repo" ] && [ -d "$tmp_repo" ] || die "mktemp gave no directory"
  git init -q --bare "$tmp_repo" || die "git init $tmp_repo failed"
  # no background gc / maintenance in a throw-away repository
  gitd "$tmp_repo" config gc.auto 0 && gitd "$tmp_repo" config maintenance.auto false || die "git config in $tmp_repo failed"
}
remote_url() {
  if [ -n "$repo_url" ]; then echo "$repo_url"; return 0; fi
  gitd "$git_dir" remote get-url origin 2>/dev/null || return 1
}
ensure_commit() {
  local commit="$1"
  # with --repo the checkout is not used (as before): the bytes come from URL
  if { [ -z "$repo_url" ] || [ "$git_dir" != "$git_dir0" ]; } && gitd "$git_dir" cat-file -e "$commit^{commit}" 2>/dev/null; then return 0; fi
  if [ "$full_fetch" = 0 ]; then
    local url err
    url="$(remote_url)" || die "$git_dir has no origin remote: use --repo URL"
    echo "fetching data commit ${commit:0:12} (trees only; blobs on demand) from $url ..." >&2
    new_tmp_repo
    gitd "$tmp_repo" remote add krhom "$url" || die "git remote add in $tmp_repo failed"
    if err="$(gitd "$tmp_repo" fetch -q --no-tags --filter=blob:none --depth 1 krhom "$commit" 2>&1)" \
       && gitd "$tmp_repo" cat-file -e "$commit^{commit}" 2>/dev/null; then
      if grep -qi 'filtering not recognized' <<<"$err"; then
        echo "fetch_library: note: the server ignored the blob filter; the whole data commit snapshot was fetched" >&2
        fetch_mode="depth-1 snapshot; the server ignored the blob filter"
      else
        fetch_mode="partial fetch of the data commit"
        [ -z "$err" ] || echo "$err" >&2
      fi
      git_dir="$tmp_repo"
      return 0
    fi
    echo "fetch_library: partial fetch refused by $url ($(head -1 <<<"$err")); falling back to a full fetch" >&2
    rm -rf "$tmp_repo"; tmp_repo=""
  fi
  if [ -n "$repo_url" ]; then new_tmp_repo; git_dir="$tmp_repo"; fetch_mode="full history of the data commit"; fi
  local remote="${repo_url:-origin}"
  echo "fetching data commit ${commit:0:12} with its history from $remote ..." >&2
  gitd "$git_dir" fetch -q --no-tags "$remote" "$commit" 2>/dev/null \
    || gitd "$git_dir" fetch -q --no-tags "$remote" "refs/heads/data-gl6:refs/remotes/krhom-fetch/data-gl6" \
    || die "cannot fetch data commit $commit from $remote"
  gitd "$git_dir" cat-file -e "$commit^{commit}" 2>/dev/null || die "data commit $commit not found in $remote (data-gl6 history)"
}

# release assets (--url): data/release_assets.tsv columns id field data_commit asset sha256 bytes
asset_bytes_total=0
fetch_url() { # fetch_url URL OUT
  if command -v curl >/dev/null 2>&1; then curl -fsSL --retry 3 --proto '=https,http,file' -o "$2" "$1"
  elif command -v wget >/dev/null 2>&1; then case "$1" in file://*) cp "${1#file://}" "$2" ;; *) wget -q -t 3 -O "$2" "$1" ;; esac
  else echo "fetch_library: need curl or wget for --url" >&2; return 1; fi
}
asset_prepare() { # asset_prepare ROW DIR: download, verify and extract the row's asset into DIR
  local row="$1" dir="$2"
  local id fld lam n prov commit costs status stl layers
  IFS=$'\t' read -r id fld lam n prov commit costs status stl layers <<<"$row"
  [ -f "$assets_tsv" ] || die "release asset list $assets_tsv not found"
  local arow; arow="$(awk -F'\t' -v id="$id" -v f="$fld" -v c="$commit" '!/^#/ && NF>=6 && $1==id && $2==f && $3==c' "$assets_tsv")"
  [ -n "$arow" ] || die "$id/$fld: no release asset for data ${commit:0:12} in $assets_tsv"
  [ "$(wc -l <<<"$arow" | tr -d ' ')" = 1 ] || die "$id/$fld: several release assets in $assets_tsv"
  local _i _f _c asset want bytes
  IFS=$'\t' read -r _i _f _c asset want bytes <<<"$arow"
  [[ "$asset" =~ ^[A-Za-z0-9._-]+$ ]] || die "$id/$fld: bad asset name '$asset' in $assets_tsv"
  echo "downloading $url_base/$asset ($bytes bytes) ..." >&2
  fetch_url "$url_base/$asset" "$dir.tar.gz" || die "cannot download $url_base/$asset"
  local got; got="$(sha_of "$dir.tar.gz")"
  [ "$got" = "$want" ] || die "$asset: sha256 $got, $assets_tsv says $want"
  asset_bytes_total=$(( asset_bytes_total + $(wc -c <"$dir.tar.gz") ))
  # listing: regular files (and their directories) at exactly the data-commit paths the row reads
  local expected listing types names extra missing e p
  expected="$( { echo "$costs"; IFS=';' read -r -a es <<<"$layers"; for e in "${es[@]}"; do p="${e%@*}"; echo "$p"; echo "$(dirname "$p")/MANIFEST.txt"; done; } | sort -u)"
  listing="$(tar -tvzf "$dir.tar.gz")" || die "$asset: not a gzip tar"
  names="$(tar -tzf "$dir.tar.gz" | sed 's#^\./##')" || die "$asset: cannot list the tar"
  types="$(awk '{print substr($1,1,1)}' <<<"$listing")"
  if grep -qv '^[-d]$' <<<"$types"; then die "$asset: tar has entries that are not regular files or directories"; fi
  extra="$(grep -v '/$' <<<"$names" | grep -vxF -f <(echo "$expected") || true)"
  [ -z "$extra" ] || die "$asset: unexpected tar entry: $(head -1 <<<"$extra")"
  if grep -q '\.\.' <<<"$names"; then die "$asset: tar entry with '..'"; fi
  missing="$(grep -vxF -f <(grep -v '/$' <<<"$names") <<<"$expected" || true)"
  [ -z "$missing" ] || die "$asset: missing $(head -1 <<<"$missing")"
  mkdir -p "$dir" || die "cannot create $dir"
  tar -C "$dir" --no-same-owner -xzf "$dir.tar.gz" || die "$asset: extraction failed"
  rm -f "$dir.tar.gz"
  echo "  $id/$fld: release asset $asset ok (sha256 ${got:0:12})" >&2
}

# src_cat COMMIT PATH: the bytes of PATH in the data commit, from git or the extracted asset
src_cat() {
  if [ -n "$url_base" ]; then cat "$staging/src/$2"; else gitd "$git_dir" show "$1:$2"; fi
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
  [ -n "$url_base" ] || ensure_commit "$commit" || die "cannot get data commit $commit"
  mkdir -p "$root" || die "cannot create the install root $root"
  staging="$(mktemp -d "$root/.staging-$id-$fld-XXXXXX")" || { staging=""; die "cannot create a staging directory under $root (not writable?)"; }
  [ -n "$staging" ] && [ -d "$staging" ] || { staging=""; die "empty or missing staging path under $root"; }
  mkdir -p "$staging/chunks" "$staging/lib/layer_a" "$staging/lib/tmpl" "$staging/lib/MANIFESTS" || die "cannot create staging subdirectories in $staging"
  [ -z "$url_base" ] || asset_prepare "$row" "$staging/src" || die "cannot prepare the release asset of $id/$fld"
  local i=0 entry path want dir manifest base got na nt man_la man_tm
  IFS=';' read -r -a entries <<<"$layers"
  for entry in "${entries[@]}"; do
    path="${entry%@*}"; want="${entry##*@}"; dir="$(dirname "$path")"
    manifest="$staging/MANIFEST.$i"
    src_cat "$commit" "$dir/MANIFEST.txt" >"$manifest" || die "no MANIFEST for $dir in $commit"
    base="$(basename "$path")"
    src_cat "$commit" "$path" >"$staging/chunks/$i" || die "cannot read $path from $commit"
    got="$(sha_of "$staging/chunks/$i")"
    [ "$got" = "$want" ] || die "$path: sha256 $got, registry says $want"
    grep -q "^$got  $base\$" "$manifest" || die "$path: sha256 $got is not the MANIFEST value"
    # tar listing: regular files and directories under layer_a/ and tmpl/ only
    local listing; listing="$(tar -tvzf "$staging/chunks/$i")" || die "$path: not a gzip tar"
    local names; names="$(tar -tzf "$staging/chunks/$i")" || die "$path: cannot list the tar"
    # the checks read here-strings / variables, never `echo | grep -q` (which fails open under
    # pipefail when grep -q exits early and echo gets SIGPIPE; review tools #3)
    local types bad_entry
    types="$(awk '{print substr($1,1,1)}' <<<"$listing")" || die "$path: cannot read the tar listing"
    if grep -qv '^[-d]$' <<<"$types"; then die "$path: tar has entries that are not regular files or directories"; fi
    bad_entry="$(grep -Ev '^(\./)?((layer_a|tmpl)/?|layer_a/a-[0-9a-f]{16}\.txt|tmpl/tmpl-[0-9a-f]{16}\.txt|\./?)$' <<<"$names" || true)"
    [ -z "$bad_entry" ] || die "$path: unexpected tar entry: $(head -1 <<<"$bad_entry")"
    na="$(grep -Ec 'layer_a/a-[0-9a-f]{16}\.txt$' <<<"$names" || true)"
    nt="$(grep -Ec 'tmpl/tmpl-[0-9a-f]{16}\.txt$' <<<"$names" || true)"
    # MANIFEST "files:" line in either word order ("layer_a 233, tmpl 471" / "249 layer_a, 233 tmpl")
    local fl; fl="$(grep -m1 '^files' "$manifest")"
    man_la="$(echo "$fl" | grep -oE '(layer_a [0-9]+|[0-9]+ layer_a)' | grep -oE '[0-9]+')"
    man_tm="$(echo "$fl" | grep -oE '(tmpl [0-9]+|[0-9]+ tmpl)' | grep -oE '[0-9]+')"
    [ "$na" = "$man_la" ] && [ "$nt" = "$man_tm" ] || die "$path: tar has $na layer_a / $nt tmpl files, MANIFEST says $man_la / $man_tm"
    grep -E '(layer_a/a-|tmpl/tmpl-)' <<<"$names" | sed 's#^\./##' | sed "s#^#$i #" >>"$staging/owners" || die "cannot write $staging/owners"
    tar -C "$staging/lib" -xzf "$staging/chunks/$i" || die "$path: extraction failed"
    cp "$manifest" "$staging/lib/MANIFESTS/$(echo "$dir" | tr '/' '_').txt" || die "cannot copy the MANIFEST of $dir"
    echo "  $id/$fld layer $i: $path ok (sha256 ${got:0:12}, $na layer_a + $nt tmpl)" >&2
    i=$((i + 1))
  done
  src_cat "$commit" "$costs" >"$staging/lib/costs.txt" || die "cannot read $costs"
  head -1 "$staging/lib/costs.txt" | grep -q '^krhom-costs v2 alg=hs-kl1' || die "$costs: wrong format"
  local status_line="complete" status_source="complete"
  if [ "$status" != "complete" ]; then
    [ -f "$repo_root/$status" ] || die "status file $status missing from this checkout"
    cp "$repo_root/$status" "$staging/lib/status.txt" || die "cannot copy $status"
    status_line="status.txt"
    status_source="$status $(sha_of "$staging/lib/status.txt")"
  fi
  # last writer of every file (registry order: later layers overwrite earlier names)
  awk '{own[$2]=$1} END {for (f in own) print own[f], f}' "$staging/owners" | sort -k2 >"$staging/final" || die "cannot write $staging/final"
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
    echo "costs_source: $costs"
    echo "status: $status_line"
    echo "status_source: $status_source"
    echo "files_layer_a: $fa"
    echo "files_tmpl: $ft"
    echo "installed: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "fetched_by: tools/fetch_library.sh @ $code_commit"
    echo "files:"
    (cd "$staging/lib" && while read -r own f; do echo "$own $(sha_of "$f") $f"; done <"$staging/final")
  } >"$staging/lib/LIBRARY.txt" || die "cannot write $staging/lib/LIBRARY.txt"
  local aside=""
  if [ -d "$dest" ]; then aside="$dest.old-$$"; mv "$dest" "$aside" || die "cannot move the old copy $dest aside"; fi
  mkdir -p "$(dirname "$dest")" || { [ -n "$aside" ] && mv "$aside" "$dest"; die "cannot create $(dirname "$dest")"; }
  if ! mv "$staging/lib" "$dest"; then
    # keep the old install: move it back before giving up
    [ -n "$aside" ] && mv "$aside" "$dest"
    die "cannot move the new copy into $dest (the old copy, if any, was kept)"
  fi
  if [ -n "$aside" ]; then rm -rf "$aside" || echo "fetch_library: warning: could not remove $aside" >&2; fi
  rm -rf "$staging" || echo "fetch_library: warning: could not remove $staging" >&2
  staging=""
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
if [ -n "$fetch_mode" ] && [ -n "$tmp_repo" ] && [ -d "$tmp_repo/objects" ]; then
  echo "downloaded $(du -sk "$tmp_repo/objects" | awk '{print $1}') KiB of git objects ($fetch_mode)" >&2
elif [ -n "$url_base" ] && [ "$asset_bytes_total" -gt 0 ]; then
  echo "downloaded $(( asset_bytes_total / 1024 )) KiB of release assets" >&2
fi
exit $rc
