# Results

Recorded outputs of the calculator. Conventions (t, q gradings, reduced/unreduced, N′) are those
of [README.md](../README.md) and [docs/CALCULATOR.md](../docs/CALCULATOR.md). The scripts, queues
and per-run logs that produced these files are on branch `archive-pre-restructure` (the former
`runs_*` directories).

## uncolored/

Checks of the uncolored engine (`engine2`, via `tools/batch_check.py`) against the KnotInfo /
LinkInfo tables in `data/`: every record with ≤ 4 strands, 1212 knots (N = 2..6) and 2138
oriented links (N = 2..5), all ok.

| file | content |
|---|---|
| `SUMMARY.md` | machine, command, meaning of the checks, coverage table |
| `knots_n2-6.tsv.gz` | all 1212 knots with ≤ 4 strands, N = 2..6 (name, crossings, strands, status, ranks) |
| `links_n2-5.tsv.gz` | the 249 links with ≤ 3 strands, N = 2..5 |
| `links4_n2-5.tsv.gz` | the 1889 four-strand links, N = 2..5 |

## colored_gl4/

∧² (color 2) homology at N′ = 2, 3, 4, computed once at box N = 4 with prune 4 from the
`l222_n4` / `l2222_n4` libraries, over the two primes F31 = F_2147483647 and F32b =
F_4294967291. Each JSON records the checks it passed (χ = decat_chi, D14 scaling, F31 = F32b,
N′ = 3 = uncolored gl_3 of the reversed braid, N′ = 2 of rank 1).

| file | content |
|---|---|
| `gl4_222.json` | three-strand knots (λ = (2,2,2)), unreduced and reduced: 62 knots |
| `gl4_222_s0.json`, `gl4_222_s1.json`, `gl4_222_s2.json` | the slice runs (16, 14 and 10 three-strand knots; `12n_767` in s1 on one prime only) |
| `gl4_222_summary.md` | one table row per three-strand knot: word, ranks at N′ = 3, 4, checks, status |
| `gl4_222_4_1_f31.json` | full `sbox` record of 4_1 over F31 (unreduced and reduced); the reference of the 4_1 anchor in `tools/randcheck.py` |
| `gl4_2222.json` | four-strand words (λ = (2,2,2,2)): Markov-stabilized and short words (m1…m5, u_s3, w41s, w52s, …) with the C2-4 / C2-7 status notes |
| `gl4_2222_c27_compare.json` | four-strand knots 6_1, 7_6, 7_7 (reduced, both primes; 7_2 pending on one prime) from the sharded `sbox` driver `tools/dev/box.sh` (session C2-7) |
