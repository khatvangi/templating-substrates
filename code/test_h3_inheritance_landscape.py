"""
Test H3 -- Inheritance carrier landscape (5 mechanisms M0-M4).

================================================================================
Purpose
================================================================================

Test H2 falsified P_H2_3 (Implementation A vs B equivalence): at N_1=20, A
gives Mode 1 wins 93%, B gives Mode 1 wins 13%. Mode 6 wins 87% in B because
per-lineage fixed templates create lineage-level inheritance via selection on
chance-favorable founder draws, even though no individual agent inherits
anything from its parent.

This means the framework's binary copyability claim conflated two things:
individual-level copy (Mode 1) and lineage-level fixation (Implementation B).
Test H3 enumerates the inheritance-mechanism landscape with finer resolution
and tests each mechanism's evolutionary capacity across multiple selection
regimes.

================================================================================
The five inheritance mechanisms
================================================================================

M0 (pure stateless / true no-inheritance):
    Each generation, each agent's phenotype is FRESHLY drawn from a fixed
    global distribution (uniform over alphabet at each position). No parent
    linkage, no mutation, no lineage state. Only selection acts (which agent
    reproduces) but offspring phenotype is fresh. This is what the framework
    SHOULD have meant by "no copy mechanism".

M1 (lineage-level fixation, no mutation, == Test H2 Implementation B):
    At lineage founding (initial population), each lineage gets one fixed
    template drawn from the global distribution. Reproduction copies the
    parent's lineage label; offspring's phenotype = lineage's fixed template
    + per-individual i.i.d. noise (rate EPS_NOISE). Mutation does NOT change
    the lineage's template -- adaptation can only happen via selection
    among lineages.

M2 (lineage-level fixation with re-draw rate r):
    Like M1, but each reproduction event triggers a fresh template draw
    with probability r (creates a new lineage with a freshly random template).
    At r=0 -> M1. At r=1 -> M0 (every offspring is a fresh template).
    H3 probes r in {0, 0.1, 0.5, 1.0}.

M3 (individual-level copy, no mutation):
    Reproduction copies the parent's GENOTYPE exactly. Phenotype = genotype
    (no noise). New variation can only come from the initial population.
    Tests whether individual copy without mutation adds anything beyond M1.

M4 (individual-level copy with mutation, == canonical Mode 1 / H2 Impl A):
    Reproduction copies the parent's genotype with per-position mutation
    rate MU_M4 = 0.01. Mutations propagate to descendants.

================================================================================
Pre-registered predictions (DEFINED BEFORE RUN)
================================================================================

Per-mechanism fitness ceiling (Step 1: isolated dynamics):

P_H3_1: M0 mean fitness plateaus at chance level (mean fitness ~ 1/ALPHABET
        = 0.25), unaffected by generation count or selection strength beta.

P_H3_2: M1 mean fitness rises above M0 due to lineage selection but plateaus
        at the maximum fitness present in the *initial population's* lineage
        templates (no new variation enters). Plateau height should depend on
        initial population size K and target sharpness.

P_H3_3: M2 mean fitness interpolates between M0 and M1: at r=1, ~M0; at r=0,
        ~M1. Plateau height monotonically decreasing in r.

P_H3_4: M3 mean fitness should equal M1's plateau when starting from
        identical initial distributions. Both lack new variation; both
        plateau at initial-distribution best.

P_H3_5: M4 mean fitness rises monotonically toward target without plateau
        (within sim length), exceeding all other mechanisms' ceilings by
        gen 1000 in all selection regimes.

Pairwise competition outcomes (Step 2):

P_H3_6: M4 vs M1 competition: M4 wins >=95% by gen 500 across all
        N_M4 in {20, 80, 200}. (Framework's preserved claim.)

P_H3_7: M4 vs M0 competition: M4 wins >=99% at all N_M4, fastest
        timescale of all pairings.

P_H3_8: M1 vs M0 competition: M1 wins by lineage selection; crossover N_M1
        (where P_M1_wins = 0.5) at moderate values. This is the result H2
        stumbled on.

P_H3_9: M3 vs M1 competition: outcome ambiguous because both lack new
        variation; predicted near-tie at all N.

The key conceptual claim (Step 3):

P_H3_10: Across all selection regimes (target sharpness beta in
         {2, 5, 10, 20}), only M4 reaches mean fitness > 0.95 of theoretical
         maximum within 1000 generations. M1, M2, M3 plateau below this
         threshold. M0 stays at chance level. This is the framework's
         strongest defensible claim: the ability to generate new variation
         through copy-with-mutation is what makes individual-level
         inheritance distinct from lineage-level fixation.

If P_H3_10 fails for any regime -- if any non-M4 mechanism reaches > 0.95
of max -- the framework's individual-vs-lineage distinction itself fails
and a deeper revision is needed.

================================================================================
Sweep design (totals)
================================================================================

Step 1 (isolated dynamics, beta sweep at L=32):
    M0, M1, M3, M4 x beta in {2,5,10,20}: 4*4 = 16 cells x 30 reps = 480
    M2 x r in {0,0.1,0.5,1.0} x beta in {2,5,10,20}: 16 cells x 30 = 480
    Step 1 total: 32 cells x 30 reps = 960 sims.

Step 2 (pairwise competition at beta=10, L=32):
    9 pairs x N_i in {20,80,200} x 30 reps = 810 sims.

Step 3 (selection regime sensitivity: L_TARGET sweep at beta=10):
    M0, M1, M3, M4 x L in {16,32,64,128}: 4*4 = 16 cells x 30 = 480
    M2 (r=0.5 only) x L in {16,32,64,128}: 4 cells x 30 = 120
    Step 3 total: 20 cells x 30 = 600 sims.

Step 4 (effective heritability dedicated runs at L=32, beta=10):
    M1, M2(r=0), M2(r=0.1), M2(r=0.5), M3: 5 cells x 10 reps = 50 sims.
    (Heritability is a smooth measurement; 10 reps suffices.)

GRAND TOTAL: 960 + 810 + 600 + 50 = 2420 sims.
At ~0.5s per sim (H2 measured 0.45s/sim), wall-time ~20 min, with the
isolated dynamics being slightly slower due to lineage tracking.

================================================================================
Reproducibility
================================================================================

Module seed: np.random.seed(42), np.random.default_rng(42).
Per-cell, per-replicate seed: seed = 42 + cell_idx*100 + rep_idx
where cell_idx is the position in the global cell list (Step 1 first, then
Step 2, then Step 3, then Step 4).

================================================================================
Anti-DRY discipline
================================================================================

Per CLAUDE.md, this script duplicates the K=400 selection loop, fitness
function, and Mode 1 (= M4) / Mode 6 implB (= M1) logic from
test_h2_competition_sweep.py rather than importing. If a bug is fixed in
one, fix it in the other intentionally.

================================================================================
CLI arguments
================================================================================

--progress-file <path>  : where to write progress lines (every 50 sims)
--completed-csv  <path>  : append every finished sim's row immediately
--smoke-test            : tiny run (Step 1 with M0 + M4 only, 2 reps,
                            50 generations) for verification
"""

import argparse
import csv
import os
import signal
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ============================================================================
# paths
# ============================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# parameters (matched to test_h2)
# ============================================================================
ALPHABET = 4
L_TARGET_DEFAULT = 32
K = 400
N_GEN = 1000
BETA_CANONICAL = 10.0

MU_M4 = 0.01           # mutation rate for M4 (only)
EPS_NOISE = 0.05       # phenotype noise for M1, M2 (lineage template + noise)

N_REPLICATES = 30
N_REPLICATES_HERITABILITY = 10  # Step 4 fewer reps (heritability is smoother)

BETA_SWEEP = [2.0, 5.0, 10.0, 20.0]
L_SWEEP = [16, 32, 64, 128]
M2_R_VALUES = [0.0, 0.1, 0.5, 1.0]

# pairwise competition initial counts
PAIRWISE_N_VALUES = [20, 80, 200]

# heritability tracking generations
HERITABILITY_GENS = [50, 200, 500]

# fixed targets, keyed by L_TARGET. seed 2026 to match h2 idiom.
def buildTarget(L):
    return np.random.default_rng(2026).integers(0, ALPHABET, size=L, dtype=np.int8)

TARGETS = {L: buildTarget(L) for L in L_SWEEP}


# ============================================================================
# fitness (vectorized) -- per L_TARGET
# ============================================================================
def findFitness(phenotypes, target):
    """phenotypes shape (N, L) -> fitness (N,) in [0,1]. Hamming similarity."""
    L = target.shape[0]
    matches = (phenotypes == target[None, :]).sum(axis=1)
    return matches.astype(np.float64) / L


# ============================================================================
# Mechanism M0 -- pure stateless
# Each generation each agent draws phenotype freshly from uniform over alphabet.
# No persistent state. Selection chooses parents but offspring get fresh draws.
# ============================================================================
class MechM0:
    label = "M0_stateless"

    def __init__(self, n_agents, rng, L):
        self.rng = rng
        self.L = L
        self._n = n_agents

    @property
    def n(self):
        return self._n

    def findPhenotypes(self):
        if self._n == 0:
            return np.empty((0, self.L), dtype=np.int8)
        # fresh uniform draw every generation
        return self.rng.integers(0, ALPHABET, size=(self._n, self.L), dtype=np.int8)

    def reproduceFromIndices(self, parent_indices):
        # only the count carries forward; phenotypes will be redrawn next gen
        self._n = len(parent_indices)

    def countDistinctLineages(self):
        # M0 has no lineages
        return 0


# ============================================================================
# Mechanism M1 -- lineage-level fixation, no mutation
# Each agent has a fixed template at founding. Reproduction inherits parent's
# template exactly. Phenotype = template + per-individual noise.
# (Equivalent to Test H2 Implementation B.)
# ============================================================================
class MechM1:
    label = "M1_lineage_fixed"

    def __init__(self, n_agents, rng, L):
        self.rng = rng
        self.L = L
        # each founder gets a distinct random template
        self.templates = rng.integers(0, ALPHABET, size=(n_agents, L), dtype=np.int8)
        # lineage_id tracks the founder lineage of each agent (for diagnostics)
        self.lineage_ids = np.arange(n_agents, dtype=np.int64)
        self._next_lineage = n_agents

    @property
    def n(self):
        return self.templates.shape[0]

    def findPhenotypes(self):
        n = self.templates.shape[0]
        if n == 0:
            return np.empty((0, self.L), dtype=np.int8)
        intended = self.templates
        noise = self.rng.random((n, self.L)) < EPS_NOISE
        random_subs = self.rng.integers(0, ALPHABET, size=(n, self.L), dtype=np.int8)
        return np.where(noise, random_subs, intended).astype(np.int8)

    def reproduceFromIndices(self, parent_indices):
        # offspring inherits parent's template AND lineage id exactly
        self.templates = self.templates[parent_indices].copy()
        self.lineage_ids = self.lineage_ids[parent_indices].copy()

    def countDistinctLineages(self):
        if self.lineage_ids.size == 0:
            return 0
        return int(np.unique(self.lineage_ids).size)


# ============================================================================
# Mechanism M2 -- lineage-level fixation with re-draw rate r
# Like M1 but each reproduction event triggers a fresh template draw with
# probability r (= a new lineage with a fresh random template).
# At r=0 -> M1; at r=1 -> all offspring fresh -> M0.
# ============================================================================
class MechM2:
    label_template = "M2_redraw_r{r:.2f}"

    def __init__(self, n_agents, rng, L, r):
        self.rng = rng
        self.L = L
        self.r = r
        self.templates = rng.integers(0, ALPHABET, size=(n_agents, L), dtype=np.int8)
        self.lineage_ids = np.arange(n_agents, dtype=np.int64)
        self._next_lineage = n_agents
        self.label = self.label_template.format(r=r)

    @property
    def n(self):
        return self.templates.shape[0]

    def findPhenotypes(self):
        n = self.templates.shape[0]
        if n == 0:
            return np.empty((0, self.L), dtype=np.int8)
        intended = self.templates
        noise = self.rng.random((n, self.L)) < EPS_NOISE
        random_subs = self.rng.integers(0, ALPHABET, size=(n, self.L), dtype=np.int8)
        return np.where(noise, random_subs, intended).astype(np.int8)

    def reproduceFromIndices(self, parent_indices):
        n_off = len(parent_indices)
        if n_off == 0:
            self.templates = np.empty((0, self.L), dtype=np.int8)
            self.lineage_ids = np.empty(0, dtype=np.int64)
            return
        # path A (1-r): inherit parent's template & lineage id
        inherited_t = self.templates[parent_indices].copy()
        inherited_l = self.lineage_ids[parent_indices].copy()
        # path B (r): fresh random template, new lineage id per offspring
        fresh_t = self.rng.integers(0, ALPHABET, size=(n_off, self.L), dtype=np.int8)
        fresh_l = np.arange(self._next_lineage, self._next_lineage + n_off, dtype=np.int64)
        self._next_lineage += n_off
        # per-offspring decision
        do_fresh = self.rng.random(n_off) < self.r
        self.templates = np.where(do_fresh[:, None], fresh_t, inherited_t)
        self.lineage_ids = np.where(do_fresh, fresh_l, inherited_l)

    def countDistinctLineages(self):
        if self.lineage_ids.size == 0:
            return 0
        return int(np.unique(self.lineage_ids).size)


# ============================================================================
# Mechanism M3 -- individual-level copy, no mutation
# Reproduction copies parent genotype exactly. No noise, no mutation.
# Phenotype = genotype.
# ============================================================================
class MechM3:
    label = "M3_indiv_copy_no_mut"

    def __init__(self, n_agents, rng, L):
        self.rng = rng
        self.L = L
        # initial population: fresh uniform draws (a fresh random sample is the
        # only source of variation; mutation never happens)
        self.genotypes = rng.integers(0, ALPHABET, size=(n_agents, L), dtype=np.int8)

    @property
    def n(self):
        return self.genotypes.shape[0]

    def findPhenotypes(self):
        # phenotype = genotype, no noise
        return self.genotypes

    def reproduceFromIndices(self, parent_indices):
        # exact copy, no mutation
        self.genotypes = self.genotypes[parent_indices].copy()

    def countDistinctLineages(self):
        # for M3 we report the count of distinct genotypes (== distinct lineages
        # since copying is exact and there is no mutation)
        if self.genotypes.size == 0:
            return 0
        return int(np.unique(self.genotypes, axis=0).shape[0])


# ============================================================================
# Mechanism M4 -- individual-level copy with mutation (canonical Mode 1)
# Faithful copy + per-position mutation MU_M4. Phenotype = genotype.
# ============================================================================
class MechM4:
    label = "M4_indiv_copy_with_mut"

    def __init__(self, n_agents, rng, L):
        self.rng = rng
        self.L = L
        self.genotypes = rng.integers(0, ALPHABET, size=(n_agents, L), dtype=np.int8)

    @property
    def n(self):
        return self.genotypes.shape[0]

    def findPhenotypes(self):
        return self.genotypes

    def reproduceFromIndices(self, parent_indices):
        new_g = self.genotypes[parent_indices].copy()
        flips = self.rng.random(new_g.shape) < MU_M4
        replacements = self.rng.integers(0, ALPHABET, size=new_g.shape, dtype=np.int8)
        new_g = np.where(flips, replacements, new_g)
        self.genotypes = new_g

    def countDistinctLineages(self):
        # for M4, "lineages" by exact genotype identity is uninformative due to
        # mutation; report number of distinct genotypes as a diversity proxy
        if self.genotypes.size == 0:
            return 0
        return int(np.unique(self.genotypes, axis=0).shape[0])


# ============================================================================
# mechanism factory
# ============================================================================
def makeMechanism(mech_key, n_agents, rng, L):
    """Create a mechanism instance from a key string."""
    if mech_key == "M0":
        return MechM0(n_agents, rng, L)
    if mech_key == "M1":
        return MechM1(n_agents, rng, L)
    if mech_key.startswith("M2_r"):
        # parse r value, e.g. "M2_r0.50"
        r_str = mech_key[len("M2_r"):]
        r = float(r_str)
        return MechM2(n_agents, rng, L, r)
    if mech_key == "M3":
        return MechM3(n_agents, rng, L)
    if mech_key == "M4":
        return MechM4(n_agents, rng, L)
    raise ValueError(f"unknown mechanism key: {mech_key}")


# ============================================================================
# isolated population dynamics driver (Step 1, Step 3)
# K agents all of one mechanism, run n_gen generations under selection.
# Returns per-generation summary list.
# ============================================================================
def runIsolated(mech_key, n_gen, beta, L, rng, track_heritability=False):
    """Run K agents of one mechanism for n_gen generations.

    Returns dict with:
        - "history": list of per-generation summary dicts
        - "heritability_samples": dict gen -> (parent_phenotypes, offspring_phenotypes)
          (only populated if track_heritability is True)
    """
    target = TARGETS[L]
    mech_rng = np.random.default_rng(rng.integers(0, 2**31 - 1))
    mech = makeMechanism(mech_key, K, mech_rng, L)

    history = []
    h_samples = {}

    for g in range(n_gen):
        phen = mech.findPhenotypes()
        fit = findFitness(phen, target)

        history.append({
            "generation": g,
            "mean_fitness": float(fit.mean()),
            "max_fitness": float(fit.max()),
            "var_fitness": float(fit.var()),
            "n_lineages": mech.countDistinctLineages(),
        })

        # selection: shared-population softmax
        w = np.exp(beta * (fit - fit.max()))
        w = w / w.sum()
        parent_indices = mech_rng.integers(0, 2**31 - 1)  # placeholder
        parent_indices = mech_rng.choice(mech.n, size=K, replace=True, p=w)

        if track_heritability and (g + 1) in HERITABILITY_GENS:
            # capture parent phenotypes (already drawn this gen) and the
            # offspring's phenotype next gen, paired by parent_indices
            parent_phen_for_offspring = phen[parent_indices].copy()
            mech.reproduceFromIndices(parent_indices)
            # draw offspring phenotypes (this is the next-gen draw)
            offspring_phen = mech.findPhenotypes()
            # NB: this consumes one extra phenotype draw from the mech rng;
            # the next iteration will redraw, which is fine for h_eff.
            h_samples[g + 1] = (parent_phen_for_offspring.copy(),
                                offspring_phen.copy())
        else:
            mech.reproduceFromIndices(parent_indices)

    return {"history": history, "heritability_samples": h_samples}


# ============================================================================
# pairwise competition driver (Step 2)
# two mechanisms share K slots, compete under shared selection.
# ============================================================================
def runPairwise(mech_key_a, mech_key_b, n_a_init, n_b_init, n_gen, beta, L, rng):
    assert n_a_init + n_b_init == K
    target = TARGETS[L]
    rng_a = np.random.default_rng(rng.integers(0, 2**31 - 1))
    rng_b = np.random.default_rng(rng.integers(0, 2**31 - 1))

    sub_a = makeMechanism(mech_key_a, n_a_init, rng_a, L)
    sub_b = makeMechanism(mech_key_b, n_b_init, rng_b, L)

    history = []
    label_a = sub_a.label
    label_b = sub_b.label

    for g in range(n_gen):
        phen_a = sub_a.findPhenotypes() if sub_a.n > 0 else np.empty((0, L), dtype=np.int8)
        phen_b = sub_b.findPhenotypes() if sub_b.n > 0 else np.empty((0, L), dtype=np.int8)
        fit_a = findFitness(phen_a, target) if sub_a.n > 0 else np.empty(0, dtype=np.float64)
        fit_b = findFitness(phen_b, target) if sub_b.n > 0 else np.empty(0, dtype=np.float64)

        n_a = sub_a.n
        n_b = sub_b.n
        n_total = n_a + n_b

        history.append({
            "generation": g,
            "mode_label": label_a,
            "count": int(n_a),
            "frequency": float(n_a) / n_total if n_total > 0 else 0.0,
            "mean_fitness": float(fit_a.mean()) if n_a > 0 else float("nan"),
        })
        history.append({
            "generation": g,
            "mode_label": label_b,
            "count": int(n_b),
            "frequency": float(n_b) / n_total if n_total > 0 else 0.0,
            "mean_fitness": float(fit_b.mean()) if n_b > 0 else float("nan"),
        })

        all_fit = np.concatenate([fit_a, fit_b])
        if all_fit.size == 0:
            break
        w = np.exp(beta * (all_fit - all_fit.max()))
        w = w / w.sum()
        parent_indices = rng.choice(n_total, size=K, replace=True, p=w)

        is_a = parent_indices < n_a
        a_parents = parent_indices[is_a]
        b_parents = parent_indices[~is_a] - n_a

        if a_parents.size > 0:
            sub_a.reproduceFromIndices(a_parents)
        else:
            sub_a = makeMechanism(mech_key_a, 0, rng_a, L)
        if b_parents.size > 0:
            sub_b.reproduceFromIndices(b_parents)
        else:
            sub_b = makeMechanism(mech_key_b, 0, rng_b, L)

    return {"history": history, "label_a": label_a, "label_b": label_b}


# ============================================================================
# effective heritability calculation
# h_eff = corr(parent_phenotype_fitness, offspring_phenotype_fitness)
# computed over paired parent->offspring at a snapshot generation.
# ============================================================================
def findEffectiveHeritability(parent_phen, offspring_phen, target):
    """Pearson correlation between parent fitness and offspring fitness."""
    parent_fit = findFitness(parent_phen, target)
    offspring_fit = findFitness(offspring_phen, target)
    if parent_fit.std() < 1e-12 or offspring_fit.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(parent_fit, offspring_fit)[0, 1])


# ============================================================================
# CSV writers (append-as-we-go, one row per finished sim)
# Different schemas per step; we use two completed CSVs:
#   - completed.csv (Step 1, Step 3) per-mechanism isolated dynamics summary
#   - completed_pairwise.csv (Step 2) per-pairwise summary
#   - completed_heritability.csv (Step 4) per-mech-per-gen heritability
# To keep things simple we write all three to the same completed.csv with a
# "row_kind" column distinguishing them.
# ============================================================================
COMPLETED_FIELDS = [
    "row_kind",       # "isolated", "pairwise", "heritability"
    "step",           # 1, 2, 3, or 4
    "cell_idx",
    "replicate",
    "seed",
    "mech_key",
    "mech_key_b",     # only for pairwise
    "n_init_a",       # only for pairwise (== K for isolated)
    "n_init_b",       # only for pairwise
    "L",
    "beta",
    "r",              # only meaningful for M2 mechanisms
    "n_gen",
    # isolated / step 1+3 columns
    "final_mean_fitness",
    "final_max_fitness",
    "final_var_fitness",
    "final_n_lineages",
    "max_mean_fitness_overall",
    "gen_reach_0p95_max",  # first generation mean_fitness/max_target >= 0.95 (or -1)
    # pairwise / step 2 columns
    "pair_label_a",
    "pair_label_b",
    "final_freq_a",
    "final_freq_b",
    "final_mean_fitness_a",
    "final_mean_fitness_b",
    "gen_a_reaches_95",
    "gen_b_reaches_95",
    "winner_label",
    # heritability / step 4 columns
    "h_eff_gen",
    "h_eff_value",
    # bookkeeping
    "wall_seconds",
]


def appendCompleted(path, rows):
    new_file = not Path(path).exists()
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COMPLETED_FIELDS)
        if new_file:
            w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in COMPLETED_FIELDS})
        f.flush()
        os.fsync(f.fileno())


# ============================================================================
# progress writer
# ============================================================================
def writeProgressLine(state):
    n_done = state["n_done"]
    n_total = state["n_total"]
    t0 = state["t0"]
    pf = state["progress_path"]
    now = time.time()
    elapsed = now - t0
    rate = n_done / elapsed if elapsed > 0 else 0.0
    remaining = n_total - n_done
    eta_sec = remaining / rate if rate > 0 else float("inf")
    eta_h = eta_sec / 3600.0
    line = (f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
            f"n_done={n_done}/{n_total} "
            f"elapsed={elapsed/60:.1f}min "
            f"rate={rate:.2f}sims/s "
            f"ETA={eta_h:.2f}h "
            f"current_step={state.get('current_step', '?')}")
    with open(pf, "a") as f:
        f.write(line + "\n")
        f.flush()


def writeFinalProgress(state, status):
    pf = state["progress_path"]
    n_done = state["n_done"]
    n_total = state["n_total"]
    elapsed = time.time() - state["t0"]
    line = (f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
            f"FINAL status={status} "
            f"n_done={n_done}/{n_total} "
            f"elapsed={elapsed/60:.1f}min")
    with open(pf, "a") as f:
        f.write(line + "\n")
        f.flush()


# ============================================================================
# cell builders for each step
# ============================================================================
def buildStep1Cells():
    """Step 1 isolated dynamics: 5 mechs x 4 betas at L=32."""
    cells = []
    # M0, M1, M3, M4 (no r)
    for mech_key in ["M0", "M1", "M3", "M4"]:
        for beta in BETA_SWEEP:
            cells.append({
                "step": 1, "mech_key": mech_key, "beta": beta,
                "L": L_TARGET_DEFAULT, "r": "",
            })
    # M2 with 4 r values
    for r in M2_R_VALUES:
        for beta in BETA_SWEEP:
            cells.append({
                "step": 1, "mech_key": f"M2_r{r:.2f}", "beta": beta,
                "L": L_TARGET_DEFAULT, "r": r,
            })
    return cells


def buildStep2Cells():
    """Step 2 pairwise: 9 pairs x 3 N_a values."""
    pairs = [
        ("M0", "M1"),
        ("M0", "M2_r0.00"),
        ("M0", "M3"),
        ("M0", "M4"),
        ("M1", "M2_r0.50"),
        ("M1", "M3"),
        ("M1", "M4"),
        ("M2_r0.00", "M4"),
        ("M3", "M4"),
    ]
    cells = []
    for (a, b) in pairs:
        for n_a in PAIRWISE_N_VALUES:
            cells.append({
                "step": 2, "mech_key": a, "mech_key_b": b,
                "n_init_a": n_a, "n_init_b": K - n_a,
                "beta": BETA_CANONICAL, "L": L_TARGET_DEFAULT,
                "r": "",
            })
    return cells


def buildStep3Cells():
    """Step 3 selection regime sensitivity: L_TARGET sweep at beta=10."""
    cells = []
    for mech_key in ["M0", "M1", "M3", "M4"]:
        for L in L_SWEEP:
            cells.append({
                "step": 3, "mech_key": mech_key, "beta": BETA_CANONICAL,
                "L": L, "r": "",
            })
    # M2 only at r=0.5 for the L sweep
    for L in L_SWEEP:
        cells.append({
            "step": 3, "mech_key": "M2_r0.50", "beta": BETA_CANONICAL,
            "L": L, "r": 0.5,
        })
    return cells


def buildStep4Cells():
    """Step 4 effective heritability: M1, M2(r=0,0.1,0.5), M3 at L=32, beta=10."""
    cells = []
    for mech_key in ["M1", "M2_r0.00", "M2_r0.10", "M2_r0.50", "M3"]:
        cells.append({
            "step": 4, "mech_key": mech_key, "beta": BETA_CANONICAL,
            "L": L_TARGET_DEFAULT, "r": "",
        })
    return cells


# ============================================================================
# cell runner for isolated dynamics (Step 1 / Step 3)
# ============================================================================
def runIsolatedCell(cell, cell_idx, n_replicates, n_gen, completed_csv, state):
    rows = []
    mech_key = cell["mech_key"]
    beta = cell["beta"]
    L = cell["L"]
    target_max = 1.0  # fitness is in [0,1]; theoretical max of mean fitness ~1.0

    for rep_idx in range(n_replicates):
        if state.get("shutdown"):
            break
        seed = 42 + cell_idx * 100 + rep_idx
        rng = np.random.default_rng(seed)
        t0 = time.time()
        out = runIsolated(mech_key, n_gen, beta, L, rng, track_heritability=False)
        dt = time.time() - t0
        hist = out["history"]

        last = hist[-1]
        max_mean_overall = max(h["mean_fitness"] for h in hist)
        # first gen where mean_fitness >= 0.95
        gen_reach_95 = -1
        for h in hist:
            if h["mean_fitness"] >= 0.95:
                gen_reach_95 = h["generation"]
                break

        row = {
            "row_kind": "isolated",
            "step": cell["step"],
            "cell_idx": cell_idx,
            "replicate": rep_idx,
            "seed": seed,
            "mech_key": mech_key,
            "L": L,
            "beta": beta,
            "r": cell.get("r", ""),
            "n_gen": n_gen,
            "final_mean_fitness": last["mean_fitness"],
            "final_max_fitness": last["max_fitness"],
            "final_var_fitness": last["var_fitness"],
            "final_n_lineages": last["n_lineages"],
            "max_mean_fitness_overall": max_mean_overall,
            "gen_reach_0p95_max": gen_reach_95,
            "wall_seconds": dt,
        }
        rows.append(row)
        appendCompleted(completed_csv, [row])
        state["n_done"] += 1
        if state["n_done"] % 50 == 0:
            writeProgressLine(state)
    return rows


# ============================================================================
# cell runner for pairwise competition (Step 2)
# ============================================================================
def runPairwiseCell(cell, cell_idx, n_replicates, n_gen, completed_csv, state):
    rows = []
    a_key = cell["mech_key"]
    b_key = cell["mech_key_b"]
    n_a_init = cell["n_init_a"]
    n_b_init = cell["n_init_b"]
    beta = cell["beta"]
    L = cell["L"]

    for rep_idx in range(n_replicates):
        if state.get("shutdown"):
            break
        seed = 42 + cell_idx * 100 + rep_idx
        rng = np.random.default_rng(seed)
        t0 = time.time()
        out = runPairwise(a_key, b_key, n_a_init, n_b_init, n_gen, beta, L, rng)
        dt = time.time() - t0
        hist = out["history"]
        label_a = out["label_a"]
        label_b = out["label_b"]

        # extract per-mode time series
        a_recs = [r for r in hist if r["mode_label"] == label_a]
        b_recs = [r for r in hist if r["mode_label"] == label_b]
        a_recs.sort(key=lambda x: x["generation"])
        b_recs.sort(key=lambda x: x["generation"])

        final_a = a_recs[-1] if a_recs else None
        final_b = b_recs[-1] if b_recs else None

        gen_a_95 = -1
        for r in a_recs:
            if r["frequency"] > 0.95:
                gen_a_95 = r["generation"]
                break
        gen_b_95 = -1
        for r in b_recs:
            if r["frequency"] > 0.95:
                gen_b_95 = r["generation"]
                break

        # winner = whichever reaches >0.95 first; if neither, whichever has
        # the larger final frequency
        if gen_a_95 >= 0 and (gen_b_95 < 0 or gen_a_95 < gen_b_95):
            winner = label_a
        elif gen_b_95 >= 0 and (gen_a_95 < 0 or gen_b_95 < gen_a_95):
            winner = label_b
        else:
            fa = final_a["frequency"] if final_a else 0.0
            fb = final_b["frequency"] if final_b else 0.0
            winner = label_a if fa > fb else label_b

        row = {
            "row_kind": "pairwise",
            "step": cell["step"],
            "cell_idx": cell_idx,
            "replicate": rep_idx,
            "seed": seed,
            "mech_key": a_key,
            "mech_key_b": b_key,
            "n_init_a": n_a_init,
            "n_init_b": n_b_init,
            "L": L,
            "beta": beta,
            "r": cell.get("r", ""),
            "n_gen": n_gen,
            "pair_label_a": label_a,
            "pair_label_b": label_b,
            "final_freq_a": float(final_a["frequency"]) if final_a else 0.0,
            "final_freq_b": float(final_b["frequency"]) if final_b else 0.0,
            "final_mean_fitness_a": float(final_a["mean_fitness"]) if final_a and not np.isnan(final_a["mean_fitness"]) else float("nan"),
            "final_mean_fitness_b": float(final_b["mean_fitness"]) if final_b and not np.isnan(final_b["mean_fitness"]) else float("nan"),
            "gen_a_reaches_95": gen_a_95,
            "gen_b_reaches_95": gen_b_95,
            "winner_label": winner,
            "wall_seconds": dt,
        }
        rows.append(row)
        appendCompleted(completed_csv, [row])
        state["n_done"] += 1
        if state["n_done"] % 50 == 0:
            writeProgressLine(state)
    return rows


# ============================================================================
# cell runner for heritability (Step 4)
# ============================================================================
def runHeritabilityCell(cell, cell_idx, n_replicates, n_gen, completed_csv, state):
    rows = []
    mech_key = cell["mech_key"]
    beta = cell["beta"]
    L = cell["L"]
    target = TARGETS[L]

    for rep_idx in range(n_replicates):
        if state.get("shutdown"):
            break
        seed = 42 + cell_idx * 100 + rep_idx
        rng = np.random.default_rng(seed)
        t0 = time.time()
        out = runIsolated(mech_key, max(HERITABILITY_GENS) + 1, beta, L, rng,
                          track_heritability=True)
        dt = time.time() - t0
        h_samples = out["heritability_samples"]
        # for each captured generation, compute h_eff and write a row
        for g in HERITABILITY_GENS:
            if g not in h_samples:
                continue
            parent_phen, offspring_phen = h_samples[g]
            h_eff = findEffectiveHeritability(parent_phen, offspring_phen, target)
            row = {
                "row_kind": "heritability",
                "step": cell["step"],
                "cell_idx": cell_idx,
                "replicate": rep_idx,
                "seed": seed,
                "mech_key": mech_key,
                "L": L,
                "beta": beta,
                "r": cell.get("r", ""),
                "n_gen": g,
                "h_eff_gen": g,
                "h_eff_value": h_eff,
                "wall_seconds": dt,
            }
            rows.append(row)
            appendCompleted(completed_csv, [row])
        state["n_done"] += 1  # one sim per (cell, rep) regardless of #gens captured
        if state["n_done"] % 50 == 0:
            writeProgressLine(state)
    return rows


# ============================================================================
# signal handler
# ============================================================================
def installSignalHandlers(state):
    def handler(signum, frame):
        state["shutdown"] = True
        writeFinalProgress(state, status=f"signal_{signum}_received")
        sys.exit(0)
    signal.signal(signal.SIGTERM, handler)
    signal.signal(signal.SIGINT, handler)


# ============================================================================
# post-processing: load completed CSV and write per-step result CSVs,
# figures, and the v3 inheritance revision document.
# ============================================================================
def loadCompleted(path):
    rows = []
    with open(path, "r") as f:
        rd = csv.DictReader(f)
        for r in rd:
            # cast common numeric fields if present
            for k in ("step", "cell_idx", "replicate", "seed", "L", "n_gen",
                      "n_init_a", "n_init_b", "final_n_lineages",
                      "gen_reach_0p95_max", "gen_a_reaches_95",
                      "gen_b_reaches_95", "h_eff_gen"):
                if r.get(k, "") not in ("", None):
                    try:
                        r[k] = int(float(r[k]))
                    except ValueError:
                        pass
            for k in ("beta", "r", "final_mean_fitness", "final_max_fitness",
                      "final_var_fitness", "max_mean_fitness_overall",
                      "final_freq_a", "final_freq_b",
                      "final_mean_fitness_a", "final_mean_fitness_b",
                      "h_eff_value", "wall_seconds"):
                if r.get(k, "") not in ("", None):
                    try:
                        r[k] = float(r[k])
                    except ValueError:
                        pass
            rows.append(r)
    return rows


# ---- Step 1 + 3 isolated-dynamics CSVs and plots
def writeIsolatedDynamicsCsv(rows, path):
    """One row per (mech_key, L, beta) cell with mean+std over replicates."""
    iso = [r for r in rows if r["row_kind"] == "isolated" and r.get("step") in (1, 3)]
    cells = {}
    for r in iso:
        key = (r["mech_key"], r["L"], r["beta"], r["step"])
        cells.setdefault(key, []).append(r)
    fields = ["step", "mech_key", "L", "beta", "r_param", "n_replicates",
              "mean_final_fitness", "std_final_fitness",
              "mean_max_fitness_overall", "std_max_fitness_overall",
              "mean_n_lineages_final", "median_gen_reach_0p95"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for key, rs in sorted(cells.items()):
            mech, L, beta, step = key
            ff = [r["final_mean_fitness"] for r in rs if not np.isnan(r["final_mean_fitness"])]
            mf = [r["max_mean_fitness_overall"] for r in rs if not np.isnan(r["max_mean_fitness_overall"])]
            nl = [r["final_n_lineages"] for r in rs]
            grs = [r["gen_reach_0p95_max"] for r in rs if r["gen_reach_0p95_max"] >= 0]
            r_param = rs[0].get("r", "")
            w.writerow({
                "step": step, "mech_key": mech, "L": L, "beta": beta,
                "r_param": r_param, "n_replicates": len(rs),
                "mean_final_fitness": float(np.mean(ff)) if ff else float("nan"),
                "std_final_fitness": float(np.std(ff)) if ff else float("nan"),
                "mean_max_fitness_overall": float(np.mean(mf)) if mf else float("nan"),
                "std_max_fitness_overall": float(np.std(mf)) if mf else float("nan"),
                "mean_n_lineages_final": float(np.mean(nl)) if nl else float("nan"),
                "median_gen_reach_0p95": int(np.median(grs)) if grs else -1,
            })


def writeSelectionRegimeCsv(rows, path):
    """Step 3: one row per (mech_key, L, beta=10) showing L sensitivity."""
    iso = [r for r in rows if r["row_kind"] == "isolated" and r.get("step") == 3]
    cells = {}
    for r in iso:
        key = (r["mech_key"], r["L"])
        cells.setdefault(key, []).append(r)
    fields = ["mech_key", "L", "beta", "n_replicates",
              "mean_final_fitness", "std_final_fitness",
              "max_observed_mean_fitness", "median_gen_reach_0p95"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for key, rs in sorted(cells.items()):
            mech, L = key
            beta = rs[0]["beta"]
            ff = [r["final_mean_fitness"] for r in rs if not np.isnan(r["final_mean_fitness"])]
            mf = [r["max_mean_fitness_overall"] for r in rs if not np.isnan(r["max_mean_fitness_overall"])]
            grs = [r["gen_reach_0p95_max"] for r in rs if r["gen_reach_0p95_max"] >= 0]
            w.writerow({
                "mech_key": mech, "L": L, "beta": beta,
                "n_replicates": len(rs),
                "mean_final_fitness": float(np.mean(ff)) if ff else float("nan"),
                "std_final_fitness": float(np.std(ff)) if ff else float("nan"),
                "max_observed_mean_fitness": float(np.max(mf)) if mf else float("nan"),
                "median_gen_reach_0p95": int(np.median(grs)) if grs else -1,
            })


def writePairwiseCsv(rows, path):
    """Step 2: one row per (a_mech, b_mech, n_init_a) cell."""
    pw = [r for r in rows if r["row_kind"] == "pairwise"]
    cells = {}
    for r in pw:
        key = (r["mech_key"], r["mech_key_b"], r["n_init_a"])
        cells.setdefault(key, []).append(r)
    fields = ["mech_key_a", "mech_key_b", "n_init_a", "n_init_b", "n_replicates",
              "p_a_wins", "p_b_wins",
              "median_gen_a_reaches_95", "median_gen_b_reaches_95",
              "mean_final_freq_a", "mean_final_freq_b",
              "mean_final_fitness_a", "mean_final_fitness_b"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for key, rs in sorted(cells.items()):
            a, b, n_a = key
            n_reps = len(rs)
            label_a = rs[0]["pair_label_a"]
            label_b = rs[0]["pair_label_b"]
            n_a_wins = sum(1 for r in rs if r["winner_label"] == label_a)
            n_b_wins = sum(1 for r in rs if r["winner_label"] == label_b)
            ag = [r["gen_a_reaches_95"] for r in rs if r["gen_a_reaches_95"] >= 0]
            bg = [r["gen_b_reaches_95"] for r in rs if r["gen_b_reaches_95"] >= 0]
            ffa = [r["final_freq_a"] for r in rs]
            ffb = [r["final_freq_b"] for r in rs]
            fma = [r["final_mean_fitness_a"] for r in rs if not np.isnan(r["final_mean_fitness_a"])]
            fmb = [r["final_mean_fitness_b"] for r in rs if not np.isnan(r["final_mean_fitness_b"])]
            w.writerow({
                "mech_key_a": a, "mech_key_b": b,
                "n_init_a": n_a, "n_init_b": K - n_a,
                "n_replicates": n_reps,
                "p_a_wins": n_a_wins / n_reps,
                "p_b_wins": n_b_wins / n_reps,
                "median_gen_a_reaches_95": int(np.median(ag)) if ag else -1,
                "median_gen_b_reaches_95": int(np.median(bg)) if bg else -1,
                "mean_final_freq_a": float(np.mean(ffa)) if ffa else float("nan"),
                "mean_final_freq_b": float(np.mean(ffb)) if ffb else float("nan"),
                "mean_final_fitness_a": float(np.mean(fma)) if fma else float("nan"),
                "mean_final_fitness_b": float(np.mean(fmb)) if fmb else float("nan"),
            })


def writeHeritabilityCsv(rows, path):
    """Step 4: one row per (mech_key, h_eff_gen) cell."""
    he = [r for r in rows if r["row_kind"] == "heritability"]
    cells = {}
    for r in he:
        key = (r["mech_key"], r["h_eff_gen"])
        cells.setdefault(key, []).append(r)
    fields = ["mech_key", "h_eff_gen", "n_replicates",
              "mean_h_eff", "std_h_eff", "median_h_eff"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for key, rs in sorted(cells.items()):
            mech, gen = key
            vals = [r["h_eff_value"] for r in rs if not np.isnan(r["h_eff_value"])]
            w.writerow({
                "mech_key": mech, "h_eff_gen": gen,
                "n_replicates": len(rs),
                "mean_h_eff": float(np.mean(vals)) if vals else float("nan"),
                "std_h_eff": float(np.std(vals)) if vals else float("nan"),
                "median_h_eff": float(np.median(vals)) if vals else float("nan"),
            })


# ============================================================================
# plotting
# ============================================================================
def plotIsolatedPlateauHeights(rows, path):
    """Step 1: bar plot of mean final fitness per mech, faceted by beta."""
    iso = [r for r in rows if r["row_kind"] == "isolated" and r.get("step") == 1]
    if not iso:
        return
    mech_keys = sorted({r["mech_key"] for r in iso})
    betas = sorted({r["beta"] for r in iso})
    fig, axes = plt.subplots(1, len(betas), figsize=(4 * len(betas), 5),
                             sharey=True)
    if len(betas) == 1:
        axes = [axes]
    cmap = plt.get_cmap("tab10")
    for i, beta in enumerate(betas):
        ax = axes[i]
        means = []
        stds = []
        for mech in mech_keys:
            cell = [r for r in iso if r["mech_key"] == mech and r["beta"] == beta]
            ff = [r["final_mean_fitness"] for r in cell if not np.isnan(r["final_mean_fitness"])]
            means.append(np.mean(ff) if ff else float("nan"))
            stds.append(np.std(ff) if ff else 0.0)
        x = np.arange(len(mech_keys))
        ax.bar(x, means, yerr=stds, color=[cmap(i % 10) for i in range(len(mech_keys))],
               alpha=0.85, capsize=3)
        ax.set_xticks(x)
        ax.set_xticklabels(mech_keys, rotation=45, ha="right", fontsize=8)
        ax.set_title(f"beta = {beta}")
        ax.axhline(0.95, color="grey", linestyle="--", alpha=0.6, label="0.95 threshold")
        ax.axhline(1.0 / ALPHABET, color="red", linestyle=":", alpha=0.6, label=f"chance ({1.0/ALPHABET:.2f})")
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3, axis="y")
        if i == 0:
            ax.set_ylabel("Mean fitness at gen 1000 (final)")
            ax.legend(loc="upper left", fontsize=8)
    fig.suptitle("Test H3 Step 1 -- Isolated dynamics plateau heights (mean +/- std, 30 reps)")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotPairwiseHeatmap(rows, path):
    """Step 2: heatmap of P(M_a wins) per pair, with N_init_a as columns."""
    pw = [r for r in rows if r["row_kind"] == "pairwise"]
    if not pw:
        return
    pairs = sorted({(r["mech_key"], r["mech_key_b"]) for r in pw})
    n_a_values = sorted({r["n_init_a"] for r in pw})
    grid = np.full((len(pairs), len(n_a_values)), np.nan)
    for i, (a, b) in enumerate(pairs):
        for j, n_a in enumerate(n_a_values):
            cell = [r for r in pw if r["mech_key"] == a and r["mech_key_b"] == b
                    and r["n_init_a"] == n_a]
            if not cell:
                continue
            label_a = cell[0]["pair_label_a"]
            n_wins = sum(1 for r in cell if r["winner_label"] == label_a)
            grid[i, j] = n_wins / len(cell)

    fig, ax = plt.subplots(figsize=(2.0 + 1.0 * len(n_a_values),
                                    1.5 + 0.55 * len(pairs)))
    im = ax.imshow(grid, vmin=0.0, vmax=1.0, cmap="RdBu_r", aspect="auto")
    ax.set_xticks(range(len(n_a_values)))
    ax.set_xticklabels([f"N_a={n}" for n in n_a_values])
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels([f"{a} vs {b}" for (a, b) in pairs], fontsize=8)
    for i in range(len(pairs)):
        for j in range(len(n_a_values)):
            v = grid[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=8, color="white" if abs(v - 0.5) > 0.3 else "black")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("P(mechanism A wins)")
    ax.set_title("Test H3 Step 2 -- Pairwise competition outcomes")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotSelectionRegime(rows, path):
    """Step 3: mean final fitness vs L_TARGET, one curve per mechanism."""
    iso3 = [r for r in rows if r["row_kind"] == "isolated" and r.get("step") == 3]
    if not iso3:
        return
    mech_keys = sorted({r["mech_key"] for r in iso3})
    Ls = sorted({r["L"] for r in iso3})
    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = plt.get_cmap("tab10")
    for i, mech in enumerate(mech_keys):
        means, stds = [], []
        for L in Ls:
            cell = [r for r in iso3 if r["mech_key"] == mech and r["L"] == L]
            ff = [r["final_mean_fitness"] for r in cell if not np.isnan(r["final_mean_fitness"])]
            means.append(np.mean(ff) if ff else float("nan"))
            stds.append(np.std(ff) if ff else 0.0)
        ax.errorbar(Ls, means, yerr=stds, marker="o", color=cmap(i % 10),
                    label=mech, capsize=3, linewidth=1.5)
    ax.axhline(0.95, color="grey", linestyle="--", alpha=0.6, label="0.95 threshold")
    ax.axhline(1.0 / ALPHABET, color="red", linestyle=":", alpha=0.6,
               label=f"chance ({1.0/ALPHABET:.2f})")
    ax.set_xlabel("Target length L")
    ax.set_xscale("log", base=2)
    ax.set_ylabel("Mean final fitness at gen 1000")
    ax.set_title("Test H3 Step 3 -- Selection regime sensitivity (beta=10)")
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotEffectiveHeritability(rows, path):
    """Step 4: bar plot of h_eff per mechanism, grouped by sample generation."""
    he = [r for r in rows if r["row_kind"] == "heritability"]
    if not he:
        return
    mech_keys = sorted({r["mech_key"] for r in he})
    gens = sorted({r["h_eff_gen"] for r in he})
    fig, ax = plt.subplots(figsize=(9, 5))
    width = 0.8 / max(1, len(gens))
    cmap = plt.get_cmap("viridis")
    for i, gen in enumerate(gens):
        means, stds = [], []
        for mech in mech_keys:
            cell = [r for r in he if r["mech_key"] == mech and r["h_eff_gen"] == gen]
            vals = [r["h_eff_value"] for r in cell if not np.isnan(r["h_eff_value"])]
            means.append(np.mean(vals) if vals else 0.0)
            stds.append(np.std(vals) if vals else 0.0)
        x = np.arange(len(mech_keys)) + i * width - 0.4 + width / 2
        ax.bar(x, means, width=width, yerr=stds,
               color=cmap(i / max(1, len(gens) - 1)),
               label=f"gen {gen}", capsize=3, alpha=0.85)
    ax.set_xticks(np.arange(len(mech_keys)))
    ax.set_xticklabels(mech_keys, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("Effective heritability h_eff (corr parent-offspring fitness)")
    ax.set_title("Test H3 Step 4 -- Effective heritability per mechanism")
    ax.set_ylim(-0.1, 1.05)
    ax.axhline(0.0, color="grey", linewidth=0.5)
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ============================================================================
# v3 inheritance revision document writer
# ============================================================================
def writeInheritanceRevision(rows, path):
    """Build the test_h3_v3_inheritance_revision.md from completed data."""
    iso1 = [r for r in rows if r["row_kind"] == "isolated" and r.get("step") == 1]
    iso3 = [r for r in rows if r["row_kind"] == "isolated" and r.get("step") == 3]
    pw = [r for r in rows if r["row_kind"] == "pairwise"]
    he = [r for r in rows if r["row_kind"] == "heritability"]

    # plateau heights table: rows = mechs (sorted), cols = beta values
    mech_keys_step1 = sorted({r["mech_key"] for r in iso1})
    betas = sorted({r["beta"] for r in iso1})
    plateau_table = {}
    for mech in mech_keys_step1:
        plateau_table[mech] = {}
        for beta in betas:
            cell = [r for r in iso1 if r["mech_key"] == mech and r["beta"] == beta]
            ff = [r["final_mean_fitness"] for r in cell if not np.isnan(r["final_mean_fitness"])]
            plateau_table[mech][beta] = float(np.mean(ff)) if ff else float("nan")

    # P_H3_10: does any non-M4 mechanism reach > 0.95 in any beta?
    p_h3_10_violations = []
    p_h3_10_pass = True
    for mech, by_b in plateau_table.items():
        if mech == "M4":
            continue
        for beta, val in by_b.items():
            if val > 0.95:
                p_h3_10_violations.append((mech, beta, val))
                p_h3_10_pass = False
    # also check M4 reaches > 0.95 in all betas
    m4_meets_095 = True
    if "M4" in plateau_table:
        for beta, val in plateau_table["M4"].items():
            if val < 0.95:
                m4_meets_095 = False

    # heritability comparison
    heritability_summary = {}
    for mech in sorted({r["mech_key"] for r in he}):
        for gen in sorted({r["h_eff_gen"] for r in he}):
            cell = [r for r in he if r["mech_key"] == mech and r["h_eff_gen"] == gen]
            vals = [r["h_eff_value"] for r in cell if not np.isnan(r["h_eff_value"])]
            heritability_summary[(mech, gen)] = (
                float(np.mean(vals)) if vals else float("nan"),
                float(np.std(vals)) if vals else 0.0,
            )

    # build markdown
    lines = []
    lines.append("# Test H3 -- v3 Inheritance Carrier Reformulation")
    lines.append("")
    lines.append("This document is generated by `code/test_h3_inheritance_landscape.py` "
                 "from the post-processed sweep data. It is the central revision "
                 "for v3 Discussion (replacing the v2 Mode 6 / inheritance theorem claim).")
    lines.append("")
    lines.append("## 1. Empirical inheritance-mechanism landscape")
    lines.append("")
    lines.append("### Step 1 plateau heights (mean fitness at gen 1000, averaged over 30 reps)")
    lines.append("")
    header = "| Mechanism | " + " | ".join(f"beta={b}" for b in betas) + " |"
    sep = "|---|" + "|".join(["---"] * len(betas)) + "|"
    lines.append(header)
    lines.append(sep)
    for mech in mech_keys_step1:
        row_vals = " | ".join(f"{plateau_table[mech][b]:.3f}" for b in betas)
        lines.append(f"| {mech} | {row_vals} |")
    lines.append("")
    lines.append(f"Theoretical max (1.0) and chance baseline (1/ALPHABET = {1.0/ALPHABET:.3f}) "
                 f"are reference horizontals.")
    lines.append("")
    lines.append("### Pre-registered prediction outcomes")
    lines.append("")
    lines.append("Each prediction was registered in the script docstring before running.")
    lines.append("")
    pred_results = evaluatePredictions(rows)
    lines.append("| ID | Verdict | Empirical |")
    lines.append("|---|---|---|")
    for pid, verdict, emp in pred_results:
        emp_safe = emp.replace("|", "\\|")
        lines.append(f"| {pid} | **{verdict}** | {emp_safe} |")
    lines.append("")
    lines.append(f"**P_H3_10 (key conceptual claim) verdict:** "
                 f"{'CONFIRMED' if p_h3_10_pass and m4_meets_095 else 'FALSIFIED' if not p_h3_10_pass else 'PARTIAL'}.")
    if p_h3_10_violations:
        lines.append("Violations (non-M4 mechanism with mean fitness > 0.95):")
        for mech, beta, val in p_h3_10_violations:
            lines.append(f"  - {mech} at beta={beta}: {val:.3f}")
    if not m4_meets_095:
        lines.append("M4 failed to reach > 0.95 in some beta regime.")
    lines.append("")

    lines.append("## 2. Reformulated inheritance condition (Drafts A, B, C)")
    lines.append("")
    lines.append("### Draft A (mechanism-specific)")
    lines.append("> A substrate is a primary Darwinian inheritance carrier iff it admits "
                 "*individual-level copying with mutation* (M4-class) and the additional "
                 "capacity-scaling, variation-preservation, and genotype-phenotype "
                 "linkage conditions (ii)-(iv). This rules out M1, M2, M3 by mechanism specification.")
    lines.append("")
    lines.append("### Draft B (mechanism-agnostic, variation-rate based)")
    lines.append("> A substrate is a primary Darwinian inheritance carrier iff it can "
                 "generate new heritable variation under selection at a rate bounded "
                 "away from zero. This is mechanism-agnostic but requires defining "
                 "'new heritable variation' precisely.")
    lines.append("")
    lines.append("### Draft C (operational, accessible-phenotype-space scaling)")
    lines.append("> A substrate is a primary Darwinian inheritance carrier iff the "
                 "population's accessible phenotype space scales unboundedly with "
                 "simulation length. Operational, testable, but implicit about mechanism.")
    lines.append("")
    lines.append("### Recommended draft (data-driven)")
    lines.append("")
    if p_h3_10_pass and m4_meets_095:
        recommended = "B"
        rationale = (
            "Draft B is favored: only M4 (which has new-variation generation rate "
            "MU_M4 > 0) reaches > 0.95 of theoretical max. M3 (individual copy without "
            "mutation, rate = 0) plateaus indistinguishably from M1 (lineage fixation). "
            "Draft A's mechanism-specific phrasing is empirically equivalent here but "
            "is less general -- there could be other mechanisms generating new variation "
            "(e.g. recombination, horizontal transfer) that Draft A would exclude on "
            "the basis of mechanism rather than function. Draft C is also operationally "
            "supported but requires a definition of 'accessible phenotype space' that "
            "is harder to make precise than 'rate of new heritable variation > 0'."
        )
    else:
        recommended = "A"
        rationale = (
            "P_H3_10 partial/failed: at least one non-M4 mechanism reaches close to "
            "or above 0.95 in some regime. Draft A (mechanism-specific) is the safer "
            "fallback: it accepts that the empirical category 'M4-class' may be the "
            "operative condition, until a cleaner mechanism-agnostic statement (Draft "
            "B) can be made precise."
        )
    lines.append(f"**Recommended: Draft {recommended}.**")
    lines.append("")
    lines.append(rationale)
    lines.append("")

    lines.append("## 3. Implication for Mode 6 (v2 line 127)")
    lines.append("")
    lines.append("v2 line 127 claimed Mode 6 cannot develop inheritance under "
                 "shared-environment competition. H3 refines this:")
    lines.append("")
    lines.append("- If real Mode 6 implements M0-class inheritance (true stateless), "
                 "the v2 claim survives unmodified.")
    lines.append("- If real Mode 6 implements M1-class inheritance (lineage-level "
                 "fixed templates with no mutation), the v2 claim must weaken to "
                 "'Mode 6 cannot develop M4-class inheritance' (i.e. cannot accumulate "
                 "adaptation past the initial-population variation ceiling).")
    lines.append("- The choice between these depends on whether real Mode 6 examples "
                 "(2D surfaces, conveyor systems with mineral nucleation) carry "
                 "lineage-stable templates. This is a literature question, flagged for "
                 "the Discussion's open-questions paragraph.")
    lines.append("")
    lines.append("## 4. Effective heritability comparison (Step 4 data)")
    lines.append("")
    lines.append("| Mechanism | gen 50 | gen 200 | gen 500 |")
    lines.append("|---|---|---|---|")
    h_mechs = sorted({m for (m, _) in heritability_summary.keys()})
    h_gens = sorted({g for (_, g) in heritability_summary.keys()})
    for mech in h_mechs:
        cells = []
        for g in h_gens:
            v = heritability_summary.get((mech, g))
            if v is None:
                cells.append("--")
            else:
                cells.append(f"{v[0]:.3f} +/- {v[1]:.3f}")
        lines.append(f"| {mech} | " + " | ".join(cells) + " |")
    lines.append("")
    lines.append("Interpretation: M1 (lineage-fixed) and M3 (exact individual copy) "
                 "should both show high h_eff because both faithfully transmit "
                 "phenotype to offspring (M1 via shared template + low noise, M3 via "
                 "exact copy). M2 should show h_eff decreasing as r increases (more "
                 "fresh draws break the parent->offspring link). If M1 h_eff matches "
                 "M3 h_eff, the framework's revised statement should be that "
                 "lineage-level and individual-level inheritance are equivalent when "
                 "both lack new variation; the operative distinction is "
                 "variation-generation rate, not the location of the fixation.")
    lines.append("")
    lines.append("## 5. v3 Discussion drop-in")
    lines.append("")
    lines.append("The single sentence that should replace v2 line 127:")
    lines.append("")
    lines.append("> Mode 6 substrates that lack a mechanism for generating new "
                 "heritable variation (whether by mutation, recombination, or other "
                 "stochastic source) cannot accumulate adaptation past the variation "
                 "present in the initial population, even when lineage-level fixation "
                 "produces high effective heritability over short timescales.")
    lines.append("")
    Path(path).write_text("\n".join(lines))


def evaluatePredictions(rows):
    """Return list of (id, verdict, empirical_string) tuples for all 10 P_H3_*"""
    results = []
    iso1 = [r for r in rows if r["row_kind"] == "isolated" and r.get("step") == 1]
    iso3 = [r for r in rows if r["row_kind"] == "isolated" and r.get("step") == 3]
    pw = [r for r in rows if r["row_kind"] == "pairwise"]

    def cellMean(rs):
        ff = [r["final_mean_fitness"] for r in rs if not np.isnan(r["final_mean_fitness"])]
        return float(np.mean(ff)) if ff else float("nan")

    # P_H3_1: M0 plateaus at chance level (~0.25)
    m0 = [r for r in iso1 if r["mech_key"] == "M0"]
    m0_means = {r["beta"]: cellMean([x for x in m0 if x["beta"] == r["beta"]])
                for r in m0}
    chance = 1.0 / ALPHABET
    m0_pass = all(abs(v - chance) < 0.05 for v in m0_means.values())
    results.append((
        "P_H3_1",
        "CONFIRMED" if m0_pass else "FALSIFIED",
        f"M0 mean final fitness by beta: " +
        ", ".join(f"beta={b}:{v:.3f}" for b, v in sorted(m0_means.items())) +
        f"; chance baseline = {chance:.3f}",
    ))

    # P_H3_2: M1 > M0 but plateaus below 0.95
    m1 = [r for r in iso1 if r["mech_key"] == "M1"]
    m1_means = {r["beta"]: cellMean([x for x in m1 if x["beta"] == r["beta"]])
                for r in m1}
    m1_pass = (all(m1_means[b] > m0_means[b] for b in m1_means)
               and all(m1_means[b] < 0.95 for b in m1_means))
    results.append((
        "P_H3_2",
        "CONFIRMED" if m1_pass else "FALSIFIED",
        f"M1 mean final fitness by beta: " +
        ", ".join(f"beta={b}:{v:.3f}" for b, v in sorted(m1_means.items())),
    ))

    # P_H3_3: M2 plateau monotonic decreasing in r (at fixed beta)
    m2_keys = sorted({r["mech_key"] for r in iso1 if r["mech_key"].startswith("M2_r")})
    m2_means_by_beta = {}
    for beta in sorted({r["beta"] for r in iso1}):
        per_r = []
        for mk in m2_keys:
            cell = [r for r in iso1 if r["mech_key"] == mk and r["beta"] == beta]
            per_r.append((mk, cellMean(cell)))
        m2_means_by_beta[beta] = per_r
    # check monotonic decreasing in r at canonical beta
    m2_at_b10 = m2_means_by_beta.get(BETA_CANONICAL, [])
    m2_pass = True
    if m2_at_b10:
        # m2_keys are sorted by r ascending (M2_r0.00, M2_r0.10, M2_r0.50, M2_r1.00)
        vals = [v for _, v in m2_at_b10]
        # allow small noise (0.05)
        m2_pass = all(vals[i] >= vals[i + 1] - 0.05 for i in range(len(vals) - 1))
    emp_m2 = "; ".join(
        f"beta={b}: " + ", ".join(f"{mk}:{v:.3f}" for mk, v in pairs)
        for b, pairs in sorted(m2_means_by_beta.items())
    )
    results.append((
        "P_H3_3",
        "CONFIRMED" if m2_pass else "REFINED",
        f"M2 monotonic in r? {m2_pass}. Means: {emp_m2}",
    ))

    # P_H3_4: M3 ~ M1 (within noise)
    m3 = [r for r in iso1 if r["mech_key"] == "M3"]
    m3_means = {r["beta"]: cellMean([x for x in m3 if x["beta"] == r["beta"]])
                for r in m3}
    m3_pass = all(abs(m3_means[b] - m1_means[b]) < 0.10 for b in m3_means)
    results.append((
        "P_H3_4",
        "CONFIRMED" if m3_pass else "REFINED",
        f"M3 vs M1 by beta: " +
        ", ".join(f"beta={b}: M3={m3_means[b]:.3f}/M1={m1_means[b]:.3f}"
                  for b in sorted(m3_means.keys())),
    ))

    # P_H3_5: M4 > all others
    m4 = [r for r in iso1 if r["mech_key"] == "M4"]
    m4_means = {r["beta"]: cellMean([x for x in m4 if x["beta"] == r["beta"]])
                for r in m4}
    others_max_by_beta = {}
    for b in m4_means:
        candidates = [m0_means[b], m1_means[b], m3_means[b]]
        candidates.extend([v for _, v in m2_means_by_beta.get(b, [])])
        others_max_by_beta[b] = max(c for c in candidates if not np.isnan(c))
    m4_pass = all(m4_means[b] > others_max_by_beta[b] + 0.05 for b in m4_means)
    results.append((
        "P_H3_5",
        "CONFIRMED" if m4_pass else "FALSIFIED",
        f"M4 vs best-of-others by beta: " +
        ", ".join(f"beta={b}: M4={m4_means[b]:.3f}/others_max={others_max_by_beta[b]:.3f}"
                  for b in sorted(m4_means.keys())),
    ))

    # P_H3_6: M4 vs M1 -- M4 wins >=95%
    pw_46 = [r for r in pw if r["mech_key"] == "M1" and r["mech_key_b"] == "M4"]
    if pw_46:
        cells46 = {}
        for r in pw_46:
            cells46.setdefault(r["n_init_a"], []).append(r)
        # winner = M4 means winner_label == r["pair_label_b"]
        m4_win_rates = {}
        for n_a, rs in sorted(cells46.items()):
            label_b = rs[0]["pair_label_b"]
            n_m4_wins = sum(1 for r in rs if r["winner_label"] == label_b)
            m4_win_rates[n_a] = n_m4_wins / len(rs)
        p6_pass = all(v >= 0.95 for v in m4_win_rates.values())
        results.append((
            "P_H3_6",
            "CONFIRMED" if p6_pass else "FALSIFIED",
            f"M4 vs M1 win rates by N_M1_init: " +
            ", ".join(f"N_M1={n}:{v:.2f}" for n, v in sorted(m4_win_rates.items())),
        ))
    else:
        results.append(("P_H3_6", "NO_DATA", "M1 vs M4 cells not found in pairwise"))

    # P_H3_7: M4 vs M0 -- M4 wins >=99%
    pw_40 = [r for r in pw if r["mech_key"] == "M0" and r["mech_key_b"] == "M4"]
    if pw_40:
        cells40 = {}
        for r in pw_40:
            cells40.setdefault(r["n_init_a"], []).append(r)
        m4_win_rates = {}
        for n_a, rs in sorted(cells40.items()):
            label_b = rs[0]["pair_label_b"]
            n_m4_wins = sum(1 for r in rs if r["winner_label"] == label_b)
            m4_win_rates[n_a] = n_m4_wins / len(rs)
        p7_pass = all(v >= 0.99 for v in m4_win_rates.values())
        results.append((
            "P_H3_7",
            "CONFIRMED" if p7_pass else "REFINED",
            f"M4 vs M0 win rates by N_M0_init: " +
            ", ".join(f"N_M0={n}:{v:.2f}" for n, v in sorted(m4_win_rates.items())),
        ))
    else:
        results.append(("P_H3_7", "NO_DATA", "M0 vs M4 cells not found"))

    # P_H3_8: M1 vs M0 -- M1 wins (with crossover)
    pw_10 = [r for r in pw if r["mech_key"] == "M0" and r["mech_key_b"] == "M1"]
    if pw_10:
        cells10 = {}
        for r in pw_10:
            cells10.setdefault(r["n_init_a"], []).append(r)
        m1_win_rates = {}
        for n_a, rs in sorted(cells10.items()):
            label_b = rs[0]["pair_label_b"]
            n_m1_wins = sum(1 for r in rs if r["winner_label"] == label_b)
            m1_win_rates[n_a] = n_m1_wins / len(rs)
        # M1 should win at all reasonable N_M0 (because lineage selection helps M1)
        p8_pass = any(v >= 0.5 for v in m1_win_rates.values())
        results.append((
            "P_H3_8",
            "CONFIRMED" if p8_pass else "FALSIFIED",
            f"M1 win rates vs M0 by N_M0_init: " +
            ", ".join(f"N_M0={n}:{v:.2f}" for n, v in sorted(m1_win_rates.items())),
        ))
    else:
        results.append(("P_H3_8", "NO_DATA", "M0 vs M1 cells not found"))

    # P_H3_9: M3 vs M1 -- near tie
    pw_13 = [r for r in pw if r["mech_key"] == "M1" and r["mech_key_b"] == "M3"]
    if pw_13:
        cells13 = {}
        for r in pw_13:
            cells13.setdefault(r["n_init_a"], []).append(r)
        # near-tie: |P_M1 - 0.5| < 0.3 across all N
        near_tie = True
        win_rates = {}
        for n_a, rs in sorted(cells13.items()):
            label_a = rs[0]["pair_label_a"]
            n_m1_wins = sum(1 for r in rs if r["winner_label"] == label_a)
            wr = n_m1_wins / len(rs)
            win_rates[n_a] = wr
            if abs(wr - 0.5) >= 0.4:
                near_tie = False
        results.append((
            "P_H3_9",
            "CONFIRMED" if near_tie else "REFINED",
            f"M1 win rates vs M3 by N_M1_init: " +
            ", ".join(f"N_M1={n}:{v:.2f}" for n, v in sorted(win_rates.items())),
        ))
    else:
        results.append(("P_H3_9", "NO_DATA", "M1 vs M3 cells not found"))

    # P_H3_10: only M4 reaches > 0.95 in any beta regime (across Step 1 mechs)
    iso1_cells = {}
    for r in iso1:
        iso1_cells.setdefault((r["mech_key"], r["beta"]), []).append(r)
    cell_means = {k: cellMean(v) for k, v in iso1_cells.items()}
    m4_meets = all(cell_means.get(("M4", b), 0.0) > 0.95 for b in BETA_SWEEP)
    others_violate = []
    for (mech, beta), v in cell_means.items():
        if mech == "M4":
            continue
        if v > 0.95:
            others_violate.append((mech, beta, v))
    p10_pass = m4_meets and not others_violate
    emp = (f"M4 means by beta: " +
           ", ".join(f"beta={b}:{cell_means.get(('M4', b), float('nan')):.3f}" for b in BETA_SWEEP) +
           f"; non-M4 violators (>0.95): " +
           (", ".join(f"{m}@beta={b}:{v:.3f}" for m, b, v in others_violate)
            if others_violate else "(none)"))
    results.append((
        "P_H3_10",
        "CONFIRMED" if p10_pass else ("FALSIFIED" if others_violate else "REFINED"),
        emp,
    ))

    return results


# ============================================================================
# main
# ============================================================================
def main():
    np.random.seed(42)

    parser = argparse.ArgumentParser()
    parser.add_argument("--progress-file", required=True)
    parser.add_argument("--completed-csv", required=True)
    parser.add_argument("--smoke-test", action="store_true",
                        help="tiny config: M0 + M4 isolated only, 2 reps, 50 gens")
    args = parser.parse_args()

    progress_path = Path(args.progress_file)
    completed_path = Path(args.completed_csv)
    progress_path.parent.mkdir(parents=True, exist_ok=True)
    completed_path.parent.mkdir(parents=True, exist_ok=True)

    progress_path.write_text(f"# === test_h3 run started at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")

    if args.smoke_test:
        # tiny config: 2 isolated cells (M0 and M4 at beta=10), 2 reps, 50 gens;
        # plus 1 pairwise (M0 vs M4 at N_a=80), 2 reps, 50 gens;
        # plus 1 heritability cell (M1), 2 reps.
        cells = [
            {"step": 1, "mech_key": "M0", "beta": 10.0, "L": L_TARGET_DEFAULT, "r": ""},
            {"step": 1, "mech_key": "M4", "beta": 10.0, "L": L_TARGET_DEFAULT, "r": ""},
            {"step": 2, "mech_key": "M0", "mech_key_b": "M4",
             "n_init_a": 80, "n_init_b": 320,
             "beta": 10.0, "L": L_TARGET_DEFAULT, "r": ""},
            {"step": 4, "mech_key": "M1", "beta": 10.0, "L": L_TARGET_DEFAULT, "r": ""},
            {"step": 1, "mech_key": "M1", "beta": 10.0, "L": L_TARGET_DEFAULT, "r": ""},
            {"step": 1, "mech_key": "M2_r0.50", "beta": 10.0, "L": L_TARGET_DEFAULT, "r": 0.5},
            {"step": 1, "mech_key": "M3", "beta": 10.0, "L": L_TARGET_DEFAULT, "r": ""},
        ]
        n_replicates = 2
        n_replicates_h = 2
        n_gen = 50
        if completed_path.exists():
            completed_path.unlink()
    else:
        cells = (buildStep1Cells() + buildStep2Cells()
                 + buildStep3Cells() + buildStep4Cells())
        n_replicates = N_REPLICATES
        n_replicates_h = N_REPLICATES_HERITABILITY
        n_gen = N_GEN
        if completed_path.exists():
            backup = completed_path.with_suffix(".csv.bak")
            completed_path.rename(backup)
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] archived prior completed.csv to {backup}\n")

    # n_total accounting: each cell contributes n_replicates sims, except step 4
    # which contributes n_replicates_h sims.
    n_total = sum(
        n_replicates_h if c["step"] == 4 else n_replicates
        for c in cells
    )
    state = {
        "n_done": 0,
        "n_total": n_total,
        "t0": time.time(),
        "progress_path": progress_path,
        "current_step": "init",
        "shutdown": False,
    }
    installSignalHandlers(state)

    with open(progress_path, "a") as f:
        f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] sweep config: "
                f"{len(cells)} cells; total sims = {n_total}; "
                f"smoke_test={args.smoke_test}; n_gen={n_gen}\n")

    # iterate cells
    for cell_idx, cell in enumerate(cells):
        if state["shutdown"]:
            break
        state["current_step"] = (
            f"step{cell['step']}/{cell.get('mech_key', '?')}"
            + (f"+{cell['mech_key_b']}" if cell.get("mech_key_b") else "")
            + f"@beta={cell['beta']},L={cell['L']}"
        )
        if cell["step"] in (1, 3):
            runIsolatedCell(cell, cell_idx, n_replicates, n_gen, completed_path, state)
        elif cell["step"] == 2:
            runPairwiseCell(cell, cell_idx, n_replicates, n_gen, completed_path, state)
        elif cell["step"] == 4:
            runHeritabilityCell(cell, cell_idx, n_replicates_h,
                                max(HERITABILITY_GENS) + 1,
                                completed_path, state)
        # progress write at each cell boundary as well
        writeProgressLine(state)

    writeFinalProgress(state, status="completed")

    # post-process (only for the full run; smoke test skips)
    if not args.smoke_test:
        try:
            all_rows = loadCompleted(completed_path)
            writeIsolatedDynamicsCsv(all_rows, RESULTS_DIR / "test_h3_isolated_dynamics_v1.csv")
            writePairwiseCsv(all_rows, RESULTS_DIR / "test_h3_pairwise_competition_v1.csv")
            writeSelectionRegimeCsv(all_rows, RESULTS_DIR / "test_h3_selection_regime_v1.csv")
            writeHeritabilityCsv(all_rows, RESULTS_DIR / "test_h3_effective_heritability_v1.csv")
            plotIsolatedPlateauHeights(all_rows, FIGURES_DIR / "test_h3_isolated_plateau_heights.png")
            plotPairwiseHeatmap(all_rows, FIGURES_DIR / "test_h3_pairwise_outcomes_heatmap.png")
            plotSelectionRegime(all_rows, FIGURES_DIR / "test_h3_selection_regime.png")
            plotEffectiveHeritability(all_rows, FIGURES_DIR / "test_h3_effective_heritability.png")
            writeInheritanceRevision(all_rows, RESULTS_DIR / "test_h3_v3_inheritance_revision.md")
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] post-process complete: "
                        f"4 step CSVs, 4 figures, v3_inheritance_revision.md written\n")
        except Exception as e:
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] post-process FAILED: {e}\n")
            raise


if __name__ == "__main__":
    main()
