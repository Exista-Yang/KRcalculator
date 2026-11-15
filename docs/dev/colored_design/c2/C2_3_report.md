# C2-3: 2-colored (2,2,2) gl_4, 3-strand knots, continuation (branch colored-2c4-k3)

Date 2026-11-13. No src/ change. Library restored from branch data-gl6 (c2/l222_n4, both primes, 233 layer-A files each, Markov corner included).
Queue: runs_c2/k3/queue2.sh (order_le10.txt = order.txt filtered to <= 10 crossings), 3-4 workers x 1 thread; one.sh now takes S from the environment.

- Knots computed: **47** (34 from C2-1 + 13 new), all 47 PASS both primes (F31 = F32b, chi_ok, scaling, N'=3 duality with reversed engine2 braid, N'=2 rank 1). 0 check failures, 0 mismatches.
- **Not computed (5 knots, 10 runs): 8_18, 10_85, 10_100, 10_148, 10_155** - NOT check failures: sbox refuses ("staged rank stage needs 10-27 GB, above the memory cap of 4500 MB"); would need a >10 GB single job. Marked as not computed (mem-cap refusal), not PASS.
- Remaining <= 10-crossing knots were still queued when this report was written (slow: 10-crossing runs take 10-60 min each).
- Cube reference, 4_1 (`krhom 4 --braid "1 -2 1 -2" --strands 3 --colors 2 --colored-method cube`, ulimit -v 10.5 GB, timeout 4 h, 1 core): **did not finish** - timed out after 14403 s with peak RSS 6.46 GB, no output. No term-by-term comparison possible (runs_c2/k3c3/).
