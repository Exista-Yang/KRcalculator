#!/usr/bin/env python3
"""Batch check of the engine2 binary against the reference tables in data/.

For every selected record (knot or oriented link) it runs the stored braid through
`engine2 STRANDS "WORD" NLIST full` (unreduced) and, with --reduced, once more with
`--reduced --load-data` on the same prepared data. Checks:

  kh      N=2 unreduced Poincare polynomial == table `kh` with q -> -q
  khred   N=2 reduced Poincare polynomial == table `kh_reduced` with q -> -q (knots only)
  chi     for every N: Euler characteristic of the unreduced homology == [N] * P(v=q^-N, z=q-q^-1)
  chired  same for reduced homology == P(v=q^-N, z=q-q^-1)
          (both sides are multiplied by z^m first, m = the largest negative power of z in P)

(P = HOMFLY-PT in the KnotInfo/Morton convention, see data/README.md; krhom puts positive
knots in negative q-degree.)  Records whose stored braid is flagged
braid_vs_invariants == "conflict" are skipped.

Output: one TSV line per record (appended, so a run can be resumed), and a summary at the end.

(engine2 is a developer binary: cargo build --release --features dev)
usage: python3 tools/batch_check.py --table data/knots_upto12.jsonl --binary target/release/engine2
         --out results.tsv [--ns 2,3,4] [--reduced] [--min-strands 2] [--max-strands 4]
         [--max-crossings 12] [--threads 1] [--timeout 1800] [--shard 0/1] [--nice 19]
"""
import argparse, json, os, re, subprocess, sys, tempfile, time
from collections import Counter, defaultdict

TERM = re.compile(r'(\d*)t\^(-?\d+)q\^(-?\d+)')


def parse_poincare(s):
    d = {}
    for c, h, q in TERM.findall(s):
        k = (int(h), int(q))
        d[k] = d.get(k, 0) + (int(c) if c else 1)
    return d


def parse_output(stdout):
    """{N: (rank, poincare dict, seconds)} and the data-preparation seconds."""
    res, prep = {}, None
    for line in stdout.splitlines():
        m = re.match(r'engine2: data ([\d.]+)s', line)
        if m:
            prep = float(m.group(1))
        m = re.match(r'N=(\d+) \S+: ([\d.]+)s rank (\d+)  (.*)$', line)
        if m:
            body = m.group(4)
            if body.strip() == '(long)':
                raise RuntimeError('engine2 truncated its output; pass "full"')
            res[int(m.group(1))] = (int(m.group(3)), parse_poincare(body), float(m.group(2)))
    return res, prep


# ---- Laurent polynomials in q as {exponent: coeff} ----
def lmul(a, b):
    out = defaultdict(int)
    for i, x in a.items():
        for j, y in b.items():
            out[i + j] += x * y
    return {k: v for k, v in out.items() if v}


def ladd(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


def qint(n):
    return {n - 1 - 2 * k: 1 for k in range(n)}


def zpow(m):
    """(q - q^-1)^m for m >= 0."""
    out = {0: 1}
    for _ in range(m):
        out = lmul(out, {1: 1, -1: -1})
    return out


def homfly_spec(homfly, n):
    """(m, z^m * P(v = q^-n, z = q - q^-1)) for terms [[a, b, c]] = sum c v^a z^b, where
    m = max(0, -min b) clears the negative powers of z that links have (z^-(c-1))."""
    m = max([0] + [-b for _, b, _ in homfly])
    out = {}
    for a, b, c in homfly:
        term = {e - n * a: c * x for e, x in zpow(b + m).items()}
        out = ladd(out, term)
    return m, out


def euler(p):
    out = defaultdict(int)
    for (h, q), k in p.items():
        out[q] += k if h % 2 == 0 else -k
    return {k: v for k, v in out.items() if v}


def table_poly(entries):
    """table [[h, q, rank]] in Bar-Natan grading -> krhom grading (h, -q)."""
    d = {}
    for h, q, k in entries:
        d[(h, -q)] = d.get((h, -q), 0) + k
    return d


def run(binary, strands, word, ns, extra, threads, timeout, nice):
    env = dict(os.environ, KRHOM_THREADS=str(threads))
    cmd = (['nice', '-n', str(nice)] if nice else []) + [binary, str(strands), word, ','.join(map(str, ns)), 'full'] + extra
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
    if p.returncode != 0:
        raise RuntimeError(f'exit {p.returncode}: {p.stderr.strip()[-300:]}')
    return parse_output(p.stdout)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--table', required=True)
    ap.add_argument('--binary', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--ns', default='2')
    ap.add_argument('--reduced', action='store_true')
    ap.add_argument('--min-strands', type=int, default=1)
    ap.add_argument('--max-strands', type=int, default=99)
    ap.add_argument('--max-crossings', type=int, default=99)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--timeout', type=float, default=1800)
    ap.add_argument('--shard', default='0/1')
    ap.add_argument('--nice', type=int, default=19)
    a = ap.parse_args()
    ns = sorted({int(x) for x in a.ns.split(',')})
    si, sm = map(int, a.shard.split('/'))

    recs = [json.loads(l) for l in open(a.table)]
    recs = [r for r in recs if a.min_strands <= r['braid']['strands'] <= a.max_strands and r['crossings'] <= a.max_crossings]
    skipped = [r['name'] for r in recs if r.get('braid_vs_invariants') == 'conflict']
    recs = [r for r in recs if r.get('braid_vs_invariants') != 'conflict']
    recs.sort(key=lambda r: (r['braid']['strands'], r['crossings'], r['name']))
    recs = [r for i, r in enumerate(recs) if i % sm == si]
    done = set()
    if os.path.exists(a.out):
        done = {l.split('\t')[0] for l in open(a.out) if l.strip()}
    print(f'{len(recs)} records selected ({len(done)} already done), skipped (braid conflict): {skipped}', flush=True)

    tmpdir = tempfile.mkdtemp(prefix='krhom_batch_')
    with open(a.out, 'a') as f:
        for r in recs:
            if r['name'] in done:
                continue
            strands = r['braid']['strands']
            word = ' '.join(str(x) for x in r['braid']['word'])
            t0 = time.time()
            fails, info = [], {}
            try:
                dpath = os.path.join(tmpdir, 'data.txt')
                extra = ['--save-data', dpath] if a.reduced else []
                un, prep = run(a.binary, strands, word, ns, extra, a.threads, a.timeout, a.nice)
                info['prep'] = prep
                info['ranks'] = {n: un[n][0] for n in ns}
                if 2 in un and r['kh'] and un[2][1] != table_poly(r['kh']):
                    fails.append('kh')
                for n in ns:
                    m, zp = homfly_spec(r['homfly'], n)
                    if lmul(zpow(m), euler(un[n][1])) != lmul(qint(n), zp):
                        fails.append(f'chi{n}')
                if a.reduced:
                    red, _ = run(a.binary, strands, word, ns, ['--reduced', '--load-data', dpath], a.threads, a.timeout, a.nice)
                    info['ranks_red'] = {n: red[n][0] for n in ns}
                    if 2 in red and 'kh_reduced' in r and red[2][1] != table_poly(r['kh_reduced']):
                        fails.append('khred')
                    for n in ns:
                        m, zp = homfly_spec(r['homfly'], n)
                        if lmul(zpow(m), euler(red[n][1])) != zp:
                            fails.append(f'chired{n}')
                status = 'ok' if not fails else 'FAIL:' + ','.join(fails)
            except subprocess.TimeoutExpired:
                status = 'timeout'
            except Exception as e:  # noqa: BLE001 - record and continue
                status = 'error'
                info['error'] = str(e)[:300]
            info['secs'] = round(time.time() - t0, 2)
            f.write('\t'.join([r['name'], str(r['crossings']), str(strands), str(r['components']), status, json.dumps(info, sort_keys=True)]) + '\n')
            f.flush()

    # summary over the whole output file
    rows = [l.rstrip('\n').split('\t') for l in open(a.out) if l.strip()]
    st = Counter(row[4].split(':')[0] for row in rows)
    print('SUMMARY', dict(st), flush=True)
    bad = [(row[0], row[4], row[5][:200]) for row in rows if row[4] != 'ok']
    for b in bad[:200]:
        print('  ', *b)
    if len(bad) > 200:
        print(f'   ... {len(bad) - 200} more')
    tim = defaultdict(list)
    for row in rows:
        try:
            tim[(int(row[2]), int(row[1]))].append(json.loads(row[5])['secs'])
        except Exception:  # noqa: BLE001
            pass
    print('timing (strands, crossings): count, median s, max s')
    for k in sorted(tim):
        v = sorted(tim[k])
        print(f'   {k}: {len(v)}, {v[len(v) // 2]:.1f}, {v[-1]:.1f}')


if __name__ == '__main__':
    main()
