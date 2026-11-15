# randcheck seed 20261115 (--inject-fault)

krhom: `/tmp/claude-0/bin/krhom`; tiers 2,3,4; jobs 2

| id | tier | word | partner | comps | N'=4 ranks u/r | checks | status | secs |
|---|---|---|---|---|---|---|---|---|
| T3-01 | 3 | `2 1 1 1 -1 -1` | markov+cancel+cancel | 1 | 6/1 | O1:PASS, O10b:WARN, O11:PASS, O12:FAIL, O2:PASS, O3:PASS, O5:FAIL, O8:PASS, O9:PASS | FAIL | 21.0 |

## Details

- T3-01 FAIL: 
  - O5 (C) FAIL: N=4 Markov/isotopy vs core u: t^0q^-4: 2 vs 1
  - O12 (C) FAIL: N=4 unknot u: t^0q^-4: 2 vs 1
  - O10b (E) WARN: N=4 negative Q coefficient -1 at t^0 q^-4
  - repro: `krhom 2,3,4 --braid 1 --strands 2 --colors 2`
  - repro: `krhom 2,3,4 --braid '2 1 1 1 -1 -1' --strands 3 --colors 2 --no-simplify --seed 201 --seed-f32b 202`
  - repro: `krhom 3 --braid '-1 -1 1 1 1 2' --strands 3`
  - repro: `krhom 2,3,4 --braid 1 --strands 2 --colors 2 --reduced`
  - repro: `krhom 2,3,4 --braid '2 1 1 1 -1 -1' --strands 3 --colors 2 --reduced --no-simplify --seed 201 --seed-f32b 202`
  - repro: `krhom 3 --braid '-1 -1 1 1 1 2' --strands 3 --reduced`

seed=20261115 n=1 PASS 0 FAIL 1 WARN 0 SKIP 0 TIMEOUT 0 ERROR 0; checks: primes 2, seeds 2, chi 4, isotopy 0, markov 1, mirror 0, duality-N3 2, reversal 0, N2 2, lewark 1, unknot 1, closed-form 6, cross-engine 0; 21s
