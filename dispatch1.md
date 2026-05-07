# Dispatch 1: Family-level Mode 3 prediction sweep (Test F)

## Working directory
`/storage/kiran-stuff/templating_framework/`

## Goal
Use the 1,232-sequence Drt3b conservation analysis to identify clade variants
defined by their secondary-gate residue signature, then predict for each
clade — under the framework's apparatus — the marginal dG fraction, period-2
autocorrelation peak, and structural information $I_\text{struct}$. Each
prediction is a falsifiable statement any structural biologist can test by
expressing the clade representative and running cDIP-seq.

This is Test F. It generalizes Test E v2 from one biological anchor to ~10.

## Background and prior results to use

Read these first (in order):
1. `results/CANONICAL_RESULTS.md` — confirms which outputs are sentinel.
2. `code/test_e_v2.py` and `results/test_e_v2_*.csv` — the working
   parameterization that gave G = 0.10 matching observed 0.1016.
3. `code/drt3b_conservation_analysis.py` and its CSV output — the residue
   identity matrix per sequence at positions E26, R253, G248, Y289, R168,
   R408, Y170, T335, T338, Y650.
4. `paper/drafts/draft_v2.md` Results section "The framework correctly
   classifies the Drt3 system" — this is what Test F extends.

## Framework-derived prediction logic

The framework's claim is structural: the universal gates (R253, G248) at
≥99% conservation enforce the cycle architecture (Mode 3, N=2); the highly
conserved primary gates (E26, Y289, Y650 at 90–94%) enforce per-state
selectivity; the variable secondary residues (R168, R408, Y170, T335, T338)
modulate fidelity ε without altering N.

For each clade variant, the framework predicts:

- $I_\text{struct} = 1 - H_2(\epsilon)$ bits per site, saturating at 1 bit
  for $\epsilon \to 0$. Any clade with intact primary gates predicts
  $I_\text{struct}$ within 0.05 bits of WT.
- Marginal dG fraction = (1/2) × P(state_A misincorporates dG) —
  determined entirely by E26 identity. Clades retaining E predict
  marginal G ≈ 0.003 (WT-like). Clades with E26 → Q substitutions
  predict ≈ 0.10 (broken state_A like Test E v2 E26Q). Clades with
  E26 → D (rare but biologically possible) predict intermediate.
- Period-2 autocorrelation peak = $(1-\epsilon)^2 + \epsilon^2$.
  Clades with intact primary gates: ≈ 0.98. Clades with degraded gates:
  ≈ 0.83 (matching Test E v2 E26Q periodicity, draft v2 line 73).
- Separation ratio against bulk-matched control: ≈ 234–1031× for intact,
  ≈ 100× for E26Q-degraded (draft v2 line 73).

The framework's strongest claim — and what makes Test F a real test — is
that secondary-residue substitutions (R168, R408, Y170, T335, T338) do NOT
shift the predicted observables outside the Mode 3 N=2 envelope. If a
secondary-residue clade variant gave (e.g.) marginal G = 0.30 or
$I_\text{struct} < 0.5$ bits, the framework's residue classification is
wrong.

## Concrete steps

### Step 1: Build clade variant table

Load `results/drt3b_conservation_analysis_residues.csv` (or whatever the
canonical output filename is — check CANONICAL_RESULTS.md). For each
sequence, the columns are residue identities at the 10 catalytic positions.

Group sequences by 10-tuple residue signature. Drop tuples with fewer than
5 representatives (keep this threshold flexible; report n per clade).
Expected outcome: 5–20 distinct clades, with the WT signature the largest.

For each clade, record:
- Clade ID (e.g. "C01_WT", "C02_E26Q-like" if any natural E→Q exists)
- Signature tuple
- Member count
- Representative organism(s)

### Step 2: Map clade signature to fidelity parameters

The Test E v2 parameterization gave WT (E26 intact) ε = 0.01 per state and
E26Q (broken state_A) p_misincorp_G = 0.45.

Define a parameterization rule (this is the framework's prediction —
document it explicitly):

- E26 = E: state_A fidelity 0.99 (WT)
- E26 = Q: state_A fidelity 0.50, dG misincorp 0.45 (Test E v2 E26Q,
  matching draft v2 reported G = 0.102 ≈ analytical 0.10)
- E26 = D: state_A fidelity 0.85, dG misincorp 0.10 (intermediate;
  carboxylate retained but geometry shifted). The draft v2 reports
  E26 conservation at 90% with conservative D substitutions accounting
  for most of the residual — so E26→D clades are EMPIRICALLY PRESENT
  in the family alignment, not synthetic. Test F's prediction for
  these clades is directly falsifiable by expressing a natural
  E26D-bearing representative.
- E26 = anything else: state_A fidelity 0.50, dG misincorp 0.45 (treat as
  E26Q-like)
- R253 = R: state_C fidelity 0.99 (WT)
- R253 ≠ R: cycle architecture broken; flag as Mode 3 disruption,
  predict $I_\text{struct} < 0.5$ bits (out-of-envelope prediction)
- G248 = G: state_C dG exclusion intact (WT)
- G248 ≠ G: predict dG fraction in C-positions rises to ~0.10
- Y289 = Y: state_C pyrimidine recognition intact
- Y289 ≠ Y: state_C selectivity degraded; predict broader marginal
- All secondary residues (R168, R408, Y170, T335, T338): no parameter
  change. The framework's prediction is that varying these does NOT shift
  observables outside the Mode 3 N=2 envelope.

Document the parameterization rule in the script header. This rule IS the
prediction.

### Step 3: Run Test E v2 simulation on each clade

For each clade, run the Test E v2 Markov-chain simulation with the clade's
parameters, sweeping L ∈ {4, 8, 16, 32, 64, 128, 256, 500} for 100 reps each.

Output per clade:
- $I_\text{struct}$ vs L
- Marginal frequencies (A, C, G, T)
- Period-2 autocorrelation peak
- Separation ratio against bulk-matched control

### Step 4: Generate prediction table

Save to `results/test_f_family_predictions_v1.csv` with columns:
- clade_id
- signature
- n_members
- predicted_marginal_G
- predicted_I_struct
- predicted_period2_peak
- predicted_separation_ratio
- envelope_classification (in-envelope / out-of-envelope)

Each row is a prediction for an experimentally testable clade
representative.

### Step 5: Make a falsifiability statement

Save `results/test_f_falsifiability_v1.md` with this structure:

> Test F predicts that for the K largest secondary-residue clades retaining
> intact primary gates (E26, R253, G248, Y289, Y650), cDIP-seq would yield
> marginal G fraction ≤ X, period-2 peak ≥ Y, $I_\text{struct} \ge Z$
> bits. Any clade observed outside this envelope falsifies the framework's
> claim that secondary residues do not affect mode classification.
> Conversely, the K' clades with naturally varying primary gates predict
> shifted observables in the directions documented in the parameterization
> rule. Confirming the predicted directions strengthens the framework;
> observing perpendicular shifts (e.g., secondary-residue clades shifting
> while primary-gate clades do not) falsifies it.

Make X, Y, Z explicit numbers, not placeholders.

## Output artifacts

- `code/test_f_family_sweep.py`
- `results/test_f_family_predictions_v1.csv`
- `results/test_f_per_clade_simulations_v1.csv` (raw simulation outputs)
- `results/test_f_falsifiability_v1.md`
- `figures/test_f_clade_predictions.png` (optional — prediction envelope plot)

Update `results/CANONICAL_RESULTS.md` to add Test F as new sentinel.

## Constraints

- Naming: all artifacts get `test_f` prefix and `_v1` suffix.
- Do NOT overwrite Test E v2 outputs.
- Document the parameterization rule (Step 2) in script comments. The rule
  is the framework's prediction; if the rule is hidden, the prediction is
  not falsifiable.
- Use 4-letter alphabet (A, C, G, T) and the Test E v2 Markov-chain
  structure exactly. Only the per-state fidelities change between clades.
- Wall-time budget: 100 reps × 8 L-values × ~15 clades × small sims should
  finish in <10 minutes. If projected runtime exceeds 15 min, dispatch
  with nohup and progress.txt.

## Reporting

Final reply should include:
1. Number of clades identified (after threshold filter)
2. Number predicted in-envelope vs out-of-envelope
3. The K-largest-clade prediction summary table (top 5 rows)
4. Path to falsifiability_v1.md
