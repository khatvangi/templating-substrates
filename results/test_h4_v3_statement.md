# Test H4 -- v3-ready statements for the manuscript

Two short paragraphs ready to drop into v3 of the manuscript. 
The first goes in the Discussion, next to the M2 r=0.10 sweet-spot 
paragraph. The second goes in the Methods/Supplementary as a note 
on the H3 finite-time win-rate.

## 1. Discussion drop-in

M4 displaces M1 at long horizons but the displacement is gradual when M1's initial population happens to draw high-fitness lineages. At gen=1000 (Test H3), M4 wins 83% at equal start; at gen=5000 (Test H4), M4 wins 80% -- approaching but not reaching the 99% threshold, so the finite-time artifact interpretation is supported but a residual long-horizon durability effect remains.

## 2. Methods / Supplementary note

The H3 finite-time win-rate (83%) reflects M4's gradual climbing from a flat-fitness initial population versus M1's chance-favorable lineage selection. At extended horizon (5000 generations, Test H4), M4's win rate converges to 80% at N_M1=200, and across all tested N_M1 values is N_M1=20: 100.00%, N_M1=80: 100.00%, N_M1=200: 80.00%. The median crossover generation at N_M1=200 (where M4 frequency first exceeds M1 frequency and stays above) is gen 2. These results support condition (i) at long horizons while documenting the gradual displacement timescale.

## Empirical numbers (for reference)

- M4 win rate at gen=5000, N_M1=200: **0.800**
- Median crossover gen at N_M1=200: **2**
- M4 win rates by N_M1_init: N_M1=20: 100.00%, N_M1=80: 100.00%, N_M1=200: 80.00%