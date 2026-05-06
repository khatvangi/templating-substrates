# Test H4 -- pre-registered prediction outcomes

Each prediction was registered in the script docstring before running. P_H4_1 is the strongest claim: M4 win rate at gen=5000 is >=99% at equal initial population (N_M1=200).

| ID | Verdict | Evidence |
|---|---|---|
| P_H4_1 | **FALSIFIED** | M4 win rate at gen=5000 with N_M1=200: 0.800 (24/30); 95% Wilson CI = [0.627, 0.905] |
| P_H4_2 | **CONFIRMED** | M1 plateau (mean fit gen 500..1000 averaged over 30 reps at N_M1=200) = 0.502 (max across reps = 0.517); reps where M1 kept climbing 500->5000 = 0/30 |
| P_H4_3 | **REFINED** | M4 final mean fitness by N_M1_init (averaged over reps where M4 still present): N_M1=20:0.973, N_M1=80:0.972, N_M1=200:0.973 |
| P_H4_4 | **CONFIRMED** | N_M1=200 crossover gen (M4 freq exceeds M1 freq and stays): median=2, 95% CI=[1, 40], reps without crossover = 2/30 |
| P_H4_5 | **CONFIRMED** | M4 win rate at N_M1=20: 1.000 (30/30) |

### P_H4_1

**Verdict:** FALSIFIED

**Evidence:** M4 win rate at gen=5000 with N_M1=200: 0.800 (24/30); 95% Wilson CI = [0.627, 0.905]

### P_H4_2

**Verdict:** CONFIRMED

**Evidence:** M1 plateau (mean fit gen 500..1000 averaged over 30 reps at N_M1=200) = 0.502 (max across reps = 0.517); reps where M1 kept climbing 500->5000 = 0/30

### P_H4_3

**Verdict:** REFINED

**Evidence:** M4 final mean fitness by N_M1_init (averaged over reps where M4 still present): N_M1=20:0.973, N_M1=80:0.972, N_M1=200:0.973

### P_H4_4

**Verdict:** CONFIRMED

**Evidence:** N_M1=200 crossover gen (M4 freq exceeds M1 freq and stays): median=2, 95% CI=[1, 40], reps without crossover = 2/30

### P_H4_5

**Verdict:** CONFIRMED

**Evidence:** M4 win rate at N_M1=20: 1.000 (30/30)
