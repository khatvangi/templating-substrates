# Test E v2 -- Drt3 Biological Anchor Classification (corrected E26Q)

Test E v2 -- corrected E26Q channel parameters from Deng et al. 2026
(P(dA)=0.80, P(dG)=0.20 at state_A in E26Q variant). Replaces the placeholder
used in v1.

## Biological context

The Drt3 system (Sharma et al. 2026, Science) is a recently-discovered
bacterial anti-phage system that produces alternating poly(GT/AC) double-
stranded DNA. The complex has two reverse transcriptases:

- Drt3a uses a conserved ACACAC region of an ncRNA as template (Mode 1).
- Drt3b synthesizes the complementary poly(AC) strand WITHOUT a nucleic acid
  template, using two amino acid residues to enforce alternation (Mode 3
  with N=2).

Sharma et al. report that the E26Q point mutation breaks one gate's
selectivity (state_A accepts dG when dGTP is available). The same protein
fold appears in AbiK, where different active-site residues yield random
DNA output -- demonstrating that "a handful of residues separates random
from sequence-specific synthesis on an identical scaffold."

## Three systems compared

  Drt3b WT      -- N=2 cyclic active site, predicted I_struct = 1 bit
  Drt3b E26Q    -- degraded N=2, one gate broken, predicted I_struct < 1 bit
  AbiK          -- same fold, no selectivity, predicted I_struct ~ 0

## Channel parameterization

  WT       state_A: P(A)=0.99, others=0.0033 each
  WT       state_C: P(C)=0.99, others=0.0033 each
  E26Q     state_A: P(A)=0.80, P(G)=0.20, P(C)=P(T)=0.0001  (Deng et al. 2026)
  E26Q     state_C: P(C)=0.99, others=0.0033 each
  AbiK            : every position uniform 0.25 each, X has no causal effect

## PASS criterion

PASS if framework apparatus correctly classifies all three systems:
  1. Drt3b WT: I_struct in [0.95, 1.0] bits, periodicity peak at lag % 2 == 0
     with value > 0.95, separation ratio > 100 vs bulk-matched control.
  2. Drt3b E26Q: I_struct < WT's value (degraded), still > 0.1 bits, separation > 5.
  3. AbiK: I_struct < 0.05 bits, separation ratio < 2.

If all three pass, the framework's biological anchor is established: the
recently-discovered Drt3 system fits exactly the Mode 3 case the framework
predicted, and the apparatus distinguishes templating from biasing-participants
in genuine biological cases.
