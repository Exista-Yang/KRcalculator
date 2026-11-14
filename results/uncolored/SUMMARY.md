# Uncolored verification (sessions U-1 and U-2) — COMPLETE for <= 4 strands

Machine: 4-core Intel Xeon 2.1GHz, 15 GB RAM, engine2 built from this branch (src/ unchanged).
Command: `tools/batch_check.py --reduced --threads 1` in 4 shards (`--shard i/4`), per-record timeout 3600 s (U-2; 1800 s in U-1),
N=2..6 (knots) and N=2..5 (links). U-2 wall time (measured with `date`): 18:03 -> 23:38 UTC, 5 h 35 min
(links <=3 strands and knots <=3 strands + first 118 four-strand knots come from U-1, ~6 min).
TSVs (gzipped, columns: name, crossings, strands, ?, status, json info with ranks):
`knots_n2-6.tsv.gz` (all 1212 knots with <=4 strands), `links_n2-5.tsv.gz` (links <=3 strands, 249),
`links4_n2-5.tsv.gz` (all 1889 four-strand links).

## Meaning of checks
- kh: N=2 unreduced Khovanov polynomial equals the table `kh` (q -> -q).
- khred: same for reduced vs `kh_reduced` (knots).
- chi / chired (N>=3, also run for N=2): Euler characteristic of unreduced / reduced homology equals [N]·P or P, with P = HOMFLY-PT at v=q^-N, z=q-q^-1.
A record is "ok" only if all of these pass for all N in the list.

## Coverage
| table | strands | N | records in table | checked | ok | fail | timeout |
|---|---|---|---|---|---|---|---|
| knots_upto12 | 2 | 2-6 | 5 | 5 | 5 | 0 | 0 |
| knots_upto12 | 3 | 2-6 | 185 | 185 | 185 | 0 | 0 |
| knots_upto12 | 4 | 2-6 | 1022 | 1022 | 1022 | 0 | 0 |
| links_upto11 | 2 | 2-5 | 6 | 6 | 6 | 0 | 0 |
| links_upto11 | 3 | 2-5 | 243 | 243 | 243 | 0 | 0 |
| links_upto11 | 4 | 2-5 | 1889 | 1889 | 1889 | 0 | 0 |

Four-strand knots by crossings: 6:1, 7:4, 8:7, 9:31, 10:74, 11:282, 12:623.
Four-strand links by crossings: 6:3, 7:8, 8:68, 9:139, 10:492, 11:1179.

Totals: knots 1212 + links 2138 = 3350 records checked, 3350 ok, 0 failures, 0 timeouts.
No records were skipped by braid conflict (none flagged among <=4 strands).
NOT covered: records with 5+ strands (excluded by --max-strands 4).
Resume: rerun the same commands; output TSVs are resumable by record name.
