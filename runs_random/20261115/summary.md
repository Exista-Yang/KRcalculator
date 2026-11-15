# randcheck seed 20261115

krhom: `/tmp/claude-0/bin/krhom`; tiers 2,3,4; jobs 2

| id | tier | word | partner | comps | N'=4 ranks u/r | checks | status | secs |
|---|---|---|---|---|---|---|---|---|
| T2-01 | 2 | `1 1 1 -1 1 -1 1 -1` |  | 2 | 36/6 | O1:PASS, O13a:PASS, O3:PASS, O4:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 2.1 |
| T2-02 | 2 | `1 1 -1 1 -1 1` |  | 2 | 36/6 | O1:PASS, O13a:PASS, O3:PASS, O4:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 1.2 |
| T2-03 | 2 | `1 1 1 -1 -1 -1` |  | 2 | 36/6 | O1:PASS, O11:PASS, O13a:PASS, O3:PASS, O4:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 1.0 |
| T2-04 | 2 | `1` |  | 1 | 6/1 | O1:PASS, O10b:PASS, O11:PASS, O13a:PASS, O3:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 0.1 |
| T2-05 | 2 | `1 1 -1` |  | 1 | 6/1 | O1:PASS, O10b:PASS, O11:PASS, O13a:PASS, O3:PASS, O4:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 0.6 |
| T2-06 | 2 | `-1 -1 -1 1 1 1 1 -1` |  | 2 | 36/6 | O1:PASS, O11:PASS, O13a:PASS, O3:PASS, O4:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 1.2 |
| T2-07 | 2 | `1 1 -1 -1 1 1 1` |  | 1 | 22/9 | O1:PASS, O10b:PASS, O11:PASS, O13a:PASS, O3:PASS, O4:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 1.8 |
| T2-08 | 2 | `1 1` |  | 2 | 36/6 | O1:PASS, O13a:PASS, O3:PASS, O4:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 0.7 |
| T3-01 | 3 | `1 1 -1 2` | markov+cancel | 1 | 6/1 | O1:PASS, O10b:PASS, O11:PASS, O12:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 2.9 |
| T3-02 | 3 | `-2 2 1 2` | markov+cancel+rotate | 1 | 6/1 | O1:PASS, O10b:PASS, O11:PASS, O12:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 25.0 |
| T3-03 | 3 | `-1 -1 2 1` | markov+rotate | 1 | 6/1 | O1:PASS, O10b:PASS, O11:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 90.7 |
| T3-04 | 3 | `-2 -1` | markov+rotate | 1 | 6/1 | O1:PASS, O10b:PASS, O11:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 2.5 |
| T3-05 | 3 | `1 -2 2 2 -1 2 -2` | markov+cancel+cancel | 2 | 36/- | O1:PASS, O11:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 161.2 |
| T3-06 | 3 | `1 -2 2 -2 1` | markov+rotate+cancel | 2 | 36/- | O1:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 61.7 |
| T3-07 | 3 | `1 2 2 1 2` | reverse | 2 | 64/- | O1:PASS, O2:PASS, O3:PASS, O7:PASS, O8:PASS, O9:PASS | PASS | 356.4 |
| T3-08 | 3 | `1 1 2 1` | reverse | 1 | 22/9 | O1:PASS, O10b:PASS, O2:PASS, O3:PASS, O7:PASS, O8:PASS, O9:PASS | PASS | 160.8 |
| T3-09 | 3 | `1 1 2 2` | mirror | 3 | 216/- | O1:PASS, O2:PASS, O3:PASS, O6:PASS, O8:PASS, O9:PASS | PASS | 39.8 |
| T3-10 | 3 | `2 -1 -1 -1` | rotate | 1 | 22/9 | O1:PASS, O10b:PASS, O2:PASS, O3:PASS, O4:PASS, O8:PASS, O9:PASS | PASS | 70.0 |
| T3-11 | 3 | `1 -2 1 -2 1 1` | reverse | 1 | 166/121 | O1:PASS, O10b:PASS, O2:PASS, O3:PASS, O7:PASS, O8:PASS, O9:PASS | PASS | 557.9 |
| T3-12 | 3 | `-1 2 2 -1 -1 -2` | reverse | 1 | 6/1 | O1:PASS, O10b:PASS, O2:PASS, O3:PASS, O7:PASS, O8:PASS, O9:PASS | PASS | 488.5 |
| T4-01 | 4 | `-1 2 -1 -1 3` | markov2+rotate | 1 | -/9 | O1:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 509.9 |
| T4-02 | 4 | `1 3 2` | markov2+rotate | 1 | 6/1 | O1:PASS, O12:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 355.6 |
| T4-03 | 4 | `-1 3 -1 -1 -2` | markov2+rotate | 1 | -/9 | O1:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 160.2 |
| T4-04 | 4 | `1 1 -3 1 -2` | markov2+rotate | 1 | -/9 | O1:PASS, O2:PASS, O3:PASS, O5:PASS, O8:PASS, O9:PASS | PASS | 177.6 |

seed=20261115 n=24 PASS 24 FAIL 0 WARN 0 SKIP 0 TIMEOUT 0 ERROR 0; checks: primes 47, seeds 31, chi 94, isotopy 10, markov 15, mirror 17, duality-N3 41, reversal 4, N2 47, lewark 22, unknot 6, closed-form 72, cross-engine 16; 1871s
