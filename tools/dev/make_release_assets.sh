#!/usr/bin/env bash
# make_release_assets.sh: write the krhom library release assets (one tarball per (library, prime))
# from the data commit, plus SHA256SUMS and the rows of data/release_assets.tsv, into a LOCAL
# directory. It uploads nothing and creates no GitHub release or tag.
#
# usage: tools/dev/make_release_assets.sh [--git-dir DIR | --repo URL] [--check] OUTDIR
#   OUTDIR      new or empty directory for the assets
#   --git-dir   git clone that has the data commit (default: this checkout); if the commit is
#               missing, a partial fetch (--filter=blob:none --depth 1) from its origin is made into
#               a temporary bare repository (the clone is not touched)
#   --repo URL  fetch from URL instead of the origin of --git-dir
#   --check     compare the written release_assets.tsv rows with data/release_assets.tsv (exit 2 on
#               a difference)
#
# Asset krhom-lib-<id>-<field>-<data7>.tar.gz holds, at their paths in the data commit, the files
# tools/fetch_library.sh --url reads for that registry row: every layer chunk, the MANIFEST.txt of
# every layer directory, and the costs file. Each chunk is checked against the registry sha256 and
# its MANIFEST before packing. The tarballs are reproducible (GNU tar: sorted names, no directory
# entries, owner 0:0, mode 0644, mtime = the data commit time; gzip -n), so --check on another
# machine with GNU tar and gzip gives the same sha256 values.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$here/../.." && pwd)"
registry="$repo_root/data/libraries.tsv"
die() { echo "make_release_assets: $*" >&2; exit 2; }
usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 1; }
git_dir="$repo_root"; repo_url=""; check=0; out=""
while [ $# -gt 0 ]; do
  case "$1" in
    --git-dir) git_dir="$2"; shift 2 ;;
    --repo) repo_url="$2"; shift 2 ;;
    --check) check=1; shift ;;
    -h|--help) usage ;;
    -*) echo "make_release_assets: unknown option $1" >&2; usage ;;
    *) [ -z "$out" ] || usage; out="$1"; shift ;;
  esac
done
[ -n "$out" ] || usage
tar --version 2>/dev/null | grep -q 'GNU tar' || die "needs GNU tar (reproducible archives)"
command -v sha256sum >/dev/null || die "needs sha256sum"
mkdir -p "$out" || die "cannot create $out"
[ -z "$(ls -A "$out")" ] || die "$out is not empty"
out="$(cd "$out" && pwd)"

tmp="$(mktemp -d)"; [ -n "$tmp" ] && [ -d "$tmp" ] || die "mktemp failed"
trap 'rm -rf "$tmp"' EXIT
src="$git_dir"
ensure_commit() {
  local commit="$1"
  if [ -z "$repo_url" ] && git -C "$src" cat-file -e "$commit^{commit}" 2>/dev/null; then return 0; fi
  if [ "$src" = "$tmp/git" ] && git -C "$src" cat-file -e "$commit^{commit}" 2>/dev/null; then return 0; fi
  local url="$repo_url"
  [ -n "$url" ] || url="$(git -C "$git_dir" remote get-url origin)" || die "no origin in $git_dir: use --repo URL"
  if [ ! -d "$tmp/git" ]; then
    git init -q --bare "$tmp/git" && git -C "$tmp/git" config gc.auto 0 && git -C "$tmp/git" remote add krhom "$url" || die "git init failed"
  fi
  src="$tmp/git"
  echo "fetching data commit ${commit:0:12} (partial) from $url ..." >&2
  git -C "$src" fetch -q --no-tags --filter=blob:none --depth 1 krhom "$commit" \
    || git -C "$src" fetch -q --no-tags krhom "$commit" || die "cannot fetch $commit from $url"
}

: >"$out/release_assets.tsv"
while IFS=$'\t' read -r id fld lam n prov commit costs status stl layers; do
  ensure_commit "$commit"
  short="${commit:0:7}"
  asset="krhom-lib-$id-$fld-$short.tar.gz"
  work="$tmp/$id-$fld"; mkdir -p "$work"
  files=("$costs")
  IFS=';' read -r -a entries <<<"$layers"
  for e in "${entries[@]}"; do
    p="${e%@*}"; want="${e##*@}"; d="$(dirname "$p")"
    for f in "$p" "$d/MANIFEST.txt"; do
      mkdir -p "$work/$(dirname "$f")"
      git -C "$src" show "$commit:$f" >"$work/$f" || die "cannot read $f from $commit"
    done
    got="$(sha256sum "$work/$p" | awk '{print $1}')"
    [ "$got" = "$want" ] || die "$p: sha256 $got, registry says $want"
    grep -q "^$got  $(basename "$p")\$" "$work/$d/MANIFEST.txt" || die "$p: not the MANIFEST sha256"
    files+=("$p" "$d/MANIFEST.txt")
  done
  mkdir -p "$work/$(dirname "$costs")"
  git -C "$src" show "$commit:$costs" >"$work/$costs" || die "cannot read $costs"
  mtime="$(git -C "$src" show -s --format=%ct "$commit")"
  printf '%s\n' "${files[@]}" | sort -u >"$tmp/list"
  tar -C "$work" --format=gnu --sort=name --no-recursion --owner=0 --group=0 --numeric-owner \
      --mode=0644 --mtime="@$mtime" -cf - -T "$tmp/list" | gzip -n -9 >"$out/$asset" || die "tar failed for $asset"
  sha="$(sha256sum "$out/$asset" | awk '{print $1}')"
  bytes="$(wc -c <"$out/$asset" | tr -d ' ')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$fld" "$commit" "$asset" "$sha" "$bytes" >>"$out/release_assets.tsv"
  echo "$asset: $bytes bytes, sha256 $sha ($(wc -l <"$tmp/list" | tr -d ' ') files)" >&2
  rm -rf "$work"
done < <(awk -F'\t' '!/^#/ && NF>=10' "$registry")
(cd "$out" && sha256sum krhom-lib-*.tar.gz >SHA256SUMS) || die "cannot write SHA256SUMS"
echo "wrote $out (SHA256SUMS, release_assets.tsv); nothing was uploaded" >&2
if [ "$check" = 1 ]; then
  if diff <(awk -F'\t' '!/^#/ && NF>=6' "$repo_root/data/release_assets.tsv") "$out/release_assets.tsv" >&2; then
    echo "release_assets.tsv: matches data/release_assets.tsv" >&2
  else
    die "the assets differ from data/release_assets.tsv (see the diff above)"
  fi
fi
