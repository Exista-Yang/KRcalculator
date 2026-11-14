# C2-S2 notes (slice2) — knots not (yet) computed

Not check failures; sbox refused or the job was stopped by the memory guard.
- 12a_850 (both primes): ran ~3 h, then sbox stopped at the final staged-rank stage: "needs about 77760 MB, above the memory cap of 12000 MB". NOT computed (mem-cap refusal), not PASS. Logs: runs_c2/k3/out_s2/12a_850_f*.err.
- 12n_466: first attempt refused by the estimate (14010 MB > 12000 MB cap); re-run with KRHOM_MEM_CAP_MB=16000, both primes PASS (real RSS ~2.5 GB).
- 12a_909: killed twice by the memory guard / cgroup OOM (RSS > 4.6 GB, 4 parallel jobs); queued last for a re-run (runs_c2/k3/out_s2/killed/).
- Memory guard (runs_c2/k3/guard_s2.sh): kills the youngest sbox when total RSS > 13.8 GB and re-queues it.
- 12n_830: F32b ran ~3 h, then sbox stopped at the final staged-rank stage: "needs about 52535 MB, above the memory cap of 16000 MB". NOT computed (mem-cap refusal); the F31 run was not attempted (same word, same stage; it would cost another ~3 h to reach the same refusal).
- 12n_708 F31: first run hit the 4 h `timeout` (exit 124, not a check failure); restarted with a 12 h timeout.
