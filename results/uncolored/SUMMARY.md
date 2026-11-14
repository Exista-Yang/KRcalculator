# Uncolored verification (session U-1) — PARTIAL

Machine: 4-core Intel Xeon 2.1GHz, 15 GB RAM. Commit of tool/engine: 07277d9 (origin/dev, includes 8fe55b5).
Command: `tools/batch_check.py --reduced --threads 1` in 4 shards (`--shard i/4`), per-record timeout 1800 s,
N=2..6 (knots) and N=2..5 (links). Wall time: about 5-6 hours (the session clock is unreliable).
TSVs: `knots_n2-6.tsv.gz`, `links_n2-5.tsv.gz` (columns: name, crossings, strands, ?, status, json info with ranks).

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
| knots_upto12 | 4 | 2-6 | 1022 | 118 (records are processed by crossing number: 1x6, 4x7, 7x8, 31x9, 72x10, 3x11 crossings; shards interleaved, so not all of each crossing number) | 118 | 0 | 0 |
| links_upto11 | 2 | 2-5 | 6 | 6 | 6 | 0 | 0 |
| links_upto11 | 3 | 2-5 | 243 | 243 | 243 | 0 | 0 |
| links_upto11 | 4 | 2-5 | 1889 | 0 | - | - | - |

Totals: 557 records checked (knots 308, links 249), 557 ok, 0 failures, 0 timeouts.
Skipped by braid conflict: the one flagged link record is a 4-strand-or-more record, not reached/selected (none among <=3 strands).
NOT covered: knots with 4 strands beyond the first 116 (mostly 10-12 crossings), all 4-strand links, and all records with 5+ strands (excluded by --max-strands 4).
4-strand records are expensive (N up to 6, 10+ crossings: minutes to >30 min each), so the full run did not fit in the session.
No failures were observed. Rerun the same commands to resume (output TSVs are resumable by record name).
