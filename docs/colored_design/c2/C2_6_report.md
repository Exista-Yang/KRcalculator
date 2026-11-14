# C2-6: 2-colored (2,2,2) gl_4, 3-strand knots, continuation (branch colored-2c4-k3)

Date 2026-11-14. No src/ change. Library restored from data-gl6 (c2/l222_n4, both primes, 233 layer-A files each, Markov corner included).
Scripts: runs_c2/k3/one2.sh (resumable, raised memory cap, ulimit per job), big.sh / big2.sh (the 5 refused knots), queue3.sh /
queue5.sh (queue: order_le10.txt first, then order.txt), memmon.sh, ckpt.sh (recompute gl4_222.json / summary, commit, push).

## Results (see the end of this file for the final counts)
Status is whatever `tools/c2/k3_c2_compare.py` printed at the last checkpoint: runs_c2/k3/compare_last.txt. Checks per knot and prime
as in C2-1 (F31 = F32b, chi_ok on 6 closures, D14 scaling, N'=3 duality with the reversed engine2 braid, N'=2 rank 1).
No check failed on any computed knot.

## The 5 refused knots
The sbox staged-rank ESTIMATE refused them (10-27 GB at cap 4.5 GB). Raising `KRHOM_MEM_CAP_MB` (12000, then 40000, then 120000: the
estimate for 10_148 was 19.2 GB, for 10_100 50.7 GB) lets them run; **all 5 computed, both primes, all checks PASS**.
Real memory (VmHWM of the sbox process, max over both primes): 8_18 1.86 GB, 10_85 2.26 GB, 10_100 3.26 GB, 10_148 1.44 GB,
10_155 1.92 GB, i.e. the estimate is 6-15x (10_100: 16x) too pessimistic, as C2-2 found. Wall (scan, 2 threads): 8_18 1323 s, 10_85 2397 s,
10_100 3443 s, 10_148 549 s, 10_155 1495 s.

## Memory behaviour of the queue
Real peak memory grows with the knot's cost in order.txt: up to ~2.8 GB per job at cost ~3e5-6e5 (10_94, 10_82, 10_106), but
10_104 (cost 7.3e5?) reached VmHWM 8.4 GB (3 threads). Three runs (10_123, 10_118, 10_99, 10_91, 10_112 attempts, kept in
runs_c2/k3/out_oom/) were killed by the cgroup OOM killer when 2-4 jobs ran in parallel (7+ GB each), or aborted at ulimit -v 6.8 GB
(10_112): these are NOT check failures, they are re-run alone with ulimit -v 12 GB (queue5.sh: cost >= 4e5 one job at a time).

## Status at this checkpoint (continuing)
- Knots computed: **61** (47 from C2-1/C2-3 + 14 new: 8_18, 10_85, 10_100, 10_148, 10_155 (the 5 refused), 10_149, 10_94, 10_159, 10_82, 10_157, 10_106, 10_99, 10_104, + 1 more listed in gl4_222.json); all **61 PASS** both primes, 0 FAIL, 0 one-prime-only.
- Largest real peaks (VmHWM): 10_99 9.8 GB (3 threads, cost 8.7e5), 10_104 8.4 GB, 10_100 3.3 GB. Jobs of cost >= 7e5 need ~8-10 GB each, so they run one at a time (ulimit -v 12 GB).
- Remaining: the rest of the <=10-crossing queue (10_116, 10_91, 10_112, 10_118, 10_109, 10_123 pending/retry) and 11-12-crossing knots;
  each heavy knot costs ~2.5-3 h (two primes sequentially), so the full 185 is far beyond one session. Runs continue; later checkpoints refresh gl4_222.json / gl4_222_summary.md.
- OOM events (not check failures): runs_c2/k3/out_oom/ (10_109, 10_123, 10_118, 10_116?/10_99 first attempt, 10_91, 10_112); out_refused_* hold the original cap refusals.
