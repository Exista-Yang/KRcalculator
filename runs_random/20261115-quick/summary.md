# randcheck seed 20261115

krhom: `/tmp/claude-0/bin/krhom`; tiers 2,3,4; jobs 2

| id | tier | word | partner | comps | N'=4 ranks u/r | checks | status | secs |
|---|---|---|---|---|---|---|---|---|
| T2-01 | 2 | `1 1 1 -1 1 -1 1 -1` |  | 2 | 36/6 | O1:PASS, O13a:PASS, O3:PASS, O4:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 2.1 |
| T2-02 | 2 | `1 1 -1 1 -1 1` |  | 2 | 36/6 | O1:PASS, O13a:PASS, O3:PASS, O4:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 1.3 |
| T2-03 | 2 | `1 1 1 -1 -1 -1` |  | 2 | 36/6 | O1:PASS, O11:PASS, O13a:PASS, O3:PASS, O4:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 1.4 |
| T2-04 | 2 | `1` |  | 1 | 6/1 | O1:PASS, O10b:PASS, O11:PASS, O13a:PASS, O3:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 0.1 |
| T3-01 | 3 | `2 1 1 1 -1 -1` | markov+cancel+cancel | 1 | 6/1 | O1:PASS, O10b:PASS, O11:PASS, O12:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 32.3 |
| T3-02 | 3 | `-1 -2` | markov+rotate+rotate | 1 | 6/1 | O1:PASS, O10b:PASS, O11:PASS, O12:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 1.6 |
| T3-03 | 3 | `-2 1 1 -1 -1 2 2` | markov+conjugate | 2 | 36/- | O1:PASS, O11:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 25.1 |
| T3-04 | 3 | `2 1 2 2` | cancel | 1 | 22/9 | O1:PASS, O10b:PASS, O2:PASS, O3:PASS, O4:PASS, O8:PASS, O9:PASS | PASS | 271.8 |
| T3-05 | 3 | `-1 -1 -2 -1` | conjugate | 1 | 22/9 | O1:PASS, O10b:PASS, O2:PASS, O3:PASS, O4:PASS, O8:PASS, O9:PASS | PASS | 192.8 |
| T3-06 | 3 | `1 1 2 2` | conjugate | 3 | 216/- | O1:PASS, O2:PASS, O3:PASS, O4:PASS, O8:PASS, O9:PASS | PASS | 85.0 |

seed=20261115 n=10 PASS 10 FAIL 0 WARN 0 SKIP 0 TIMEOUT 0 ERROR 0; checks: primes 21, seeds 13, chi 42, isotopy 6, markov 5, mirror 8, duality-N3 18, reversal 0, N2 21, lewark 10, unknot 4, closed-form 35, cross-engine 8; 484s
