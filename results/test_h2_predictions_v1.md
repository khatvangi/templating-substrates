# Test H2 -- Pre-registered prediction results (v1)

Pre-registered in `code/test_h2_competition_sweep.py` header BEFORE simulation runs.

## Predictions table

| ID      | Verdict      | Empirical observation |
|---------|--------------|------------------------|
| P_H2_1 | **FALSIFIED** | P(Mode 1 wins) sequence by N_1: 20:0.93, 40:0.90, 60:1.00, 80:0.97, 100:0.97, 120:1.00, 160:0.97, 200:1.00, 240:0.97, 280:1.00, 320:1.00, 360:1.00, 380:1.00; monotonic_within_noise=True; p@N_1=20=0.9333333333333333; p@N_1=80=0.9666666666666667 |
| P_H2_2 | **INCONCLUSIVE** | Implementation A boundary: N_1 at P_wins=0.5 is ~None, at P_wins=0.95 is ~49.999999999999986 |
| P_H2_3 | **FALSIFIED** | max per-cell \|P_A - P_B\| = 0.800 across 6 shared N_1 cells; per-cell A vs B: N_1=20:A=0.93/B=0.13, N_1=40:A=0.90/B=0.23, N_1=80:A=0.97/B=0.53, N_1=120:A=1.00/B=0.70, N_1=160:A=0.97/B=0.83, N_1=200:A=1.00/B=0.83 |
| P_H2_4 | **REFINED** | crossover N_1 at h=0.5 is ~None; P(Mode 1 wins) at N_1=200, h=1.0 is 0.9 |
| P_H2_5 | **REFINED** | crossover N_1 (where P_wins=0.5) by h: h=0.00:N_1~None, h=0.25:N_1~None, h=0.50:N_1~None, h=0.75:N_1~None, h=1.00:N_1~None |

## Notes

- Sweep totals: Implementation A = 13 N_1 cells x 30 reps; Implementation B = 6 cells x 30 reps; Implementation C = 5 h x 11 N_1 cells x 30 reps.
- K = 400, L_TARGET = 32, N_GEN = 1000, BETA = 10.0, MU_MODE1 = 0.01, EPS_NOISE_6 = 0.05
- 95% CI from Wilson score on the per-cell binomial.
- Seeding: seed = 42 + cell_idx*100 + rep_idx (unique per (cell, rep) globally).