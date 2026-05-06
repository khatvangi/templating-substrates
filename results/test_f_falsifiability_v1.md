# Test F -- falsifiability statement (v1)

## Predictions

Test F predicts that for the **K = 33** largest secondary-residue clades retaining
intact universal gates (R253, G248) and primary gates (E26, Y289, Y650),
cDIP-seq on the natural family representatives would yield:

- marginal G fraction <= **X = 0.0033** at L=64
- period-2 autocorrelation peak >= **Y = 0.9793** at L=64
- $I_\text{struct}$ >= **Z = 0.9996** bits at L=64

(Numbers above are the worst-case observable across all in-envelope clades
in the simulation; the framework predicts every member clade meets each
threshold.)

Any clade observed outside this envelope -- e.g. a secondary-residue-only
representative for which marginal G > 0.003 or $I_\text{struct}$ < 1.000
bits -- **falsifies** the framework's claim that secondary residues
(R168, R408, Y170, T335, T338) do not affect mode classification.

Conversely, the **K' = 9** clades with naturally varying primary gates
predict shifted observables in the directions documented in the
parameterization rule (see code/test_f_family_sweep.py module docstring):

- marginal G in the simulated range [0.0033, 0.0517]
- period-2 peak in the simulated range [0.8548, 0.9810]

Confirming the predicted **direction** for these clades (G shifts up,
period-2 peak shifts down for E26->Q/D substitutions) strengthens the
framework. Observing **perpendicular** shifts -- secondary-residue clades
shifting while primary-gate clades do not -- falsifies it.

## Universal-gate disrupted clades

0 clades (out of 42 total) carry universal-gate
substitutions (R253 != R or G248 != G). Framework predicts these are
out-of-Mode-3 entirely: $I_\text{struct}$ < 0.5 bits, period-2 peak ~0.25
(no cycle structure). These are the strongest pre-registered out-of-envelope
predictions, but they are not testable in the natural alignment because
both universal gates are essentially fixed across all 1,232 sequences.
Confirming the prediction would require site-directed mutagenesis
(R253A or G248A) and cDIP-seq, not family sampling.

## Y650 is primary but does not shift elongation observables

Three clades in the alignment carry a Y650F substitution as their only
primary-gate change (clades labelled `_PRIM-Y650F` in the CSV). The
parameterization rule documents that Y650 is the C-terminal protein-priming
tyrosine -- it controls *initiation* rather than the per-state elongation
channel. So the framework predicts that Y650-only-substituted clades stay
in the WT-like elongation envelope. All three Y650F clades do (marginal
G ~0.003, period-2 peak ~0.98); a Y650F clade observed *outside* the
elongation envelope would be a partial counter-example for the channel
parameterization rule even though it would not falsify the universal/primary
partition.

## Experimental test

Each row in `test_f_family_predictions_v1.csv` names a representative
sequence from the Deng et al. 2026 family alignment. Express the
representative, run cDIP-seq as in Sharma et al. 2026, and check whether
each observable falls within the row's predicted envelope. A single
violation in the predicted direction across multiple clades constitutes
a meaningful counter-example.

## Pre-registered numerical thresholds

| quantity                | predicted bound                |
|-------------------------|--------------------------------|
| marginal G (in-envel.)  | <= 0.0033                     |
| period-2 peak (in-env.) | >= 0.9793                     |
| I_struct (in-env.)      | >= 0.9996 bits                |
| marginal G (primary)    | in [0.0033, 0.0517] |
| period-2 peak (primary) | in [0.8548, 0.9810] |
