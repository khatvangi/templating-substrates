## Test D v2 -- corrected Mode 6 (no template inheritance)

# Test D -- Generational information / population dynamics

## Goal

Test the framework's deepest claim (P6): only Modes 1 and 2 can serve as
substrates for inherited information across generations under selection.
Modes 3, 4, 5, 6 hit information-capacity ceilings.

## Setup

- Population size: 200 agents
- Generations: 1000
- Selection: probability proportional to exp(10 * fitness)
- Fitness: fraction of phenotype positions matching a fixed random target of
  length 32 over 4-letter alphabet
- Five modes simulated, each with mode-specific inheritance and mutation rules
- Five 2nd-order-matched bulk controls (per Gemini's review)

## Mode-specific rules

  Mode 1: 1D template, faithful copy with point mutation rate 0.01
  Mode 3 (N=2): two-state cyclic active site, mutation flips one state's preference
  Mode 4: single conformer (1 nucleotide), all positions are the same
  Mode 5 (N_modules=8): 8 modules; positions 0-7 templated, 8-31 untemplated random
  Mode 6: 2D surface templating with lossy inheritance (10% per-position copy error)

## 2nd-order bulk controls

For each mode, a parallel "bulk" version where phenotypes are drawn from a
2nd-order Markov chain matching the mode's empirical 1st- and 2nd-order
statistics, but DECOUPLED from individual parental templates. Bulk controls
should not climb the fitness landscape; their fitness reflects only what
2nd-order statistics provide.

## Predicted ceilings (mean fitness at generation 1000)

  Mode 1:        ~0.95+ (climbs toward perfect match)
  Mode 3 (N=2):  ~0.40 (best alternating pattern matches ~half of random target)
  Mode 4:        ~0.30 (best single nucleotide matches ~1/4 + small selection lift)
  Mode 5 (N=8):  ~0.4375 (8 perfect templated + 24 random at 0.25 = 14/32)
  Mode 6:        ~0.40 (lossy inheritance erodes accumulated information)
  Bulk controls: ~0.30 (chance + small lift from 2nd-order matching)

## PASS criterion

PASS if:
  1. Mode 1 mean fitness > 0.85 by gen 1000
  2. Mode 3 plateau <= 0.40
  3. Mode 4 plateau <= 0.30
  4. Mode 5 plateau in [N_modules/L_target - 0.10, N_modules/L_target + 0.10]
  5. Mode 6 plateau <= 0.40
  6. All bulk controls plateau <= 0.30
  7. Ordering: Mode_1 > Mode_5 > Mode_6 ~ Mode_3 > Mode_4 at gen 1000

If all 7 hold, the framework's deepest claim is empirically validated:
inheritance-capable templating (Mode 1) accumulates adaptation while
information-bounded templating modes plateau at their capacities.

## Why v2 exists: v1 modeled Mode 6 as lossy 1D copy with mu=0.10, which gave it Mode-1-like behavior (fitness 0.634). The framework's actual claim is that 2D surfaces lack a complementarity-based copy mechanism. v2 implements this honestly: each offspring surface is randomly initialized, blocking generational information transfer. Compare v1 vs v2 results to see how Mode 6 specification matters.
