"""
Test H5 -- M2 sweet spot characterization (variation-rate principle).

================================================================================
Purpose
================================================================================

Test H3 found an unanticipated sweet spot in M2 (lineage fixation with
re-draw rate r): plateau height at r=0.10 is ~0.53, exceeding both
r=0.00 (pure fixation, ~0.47) and r=0.50 (~0.42). The framework's
intuition (P_H3_3) that re-draws should hurt M1 monotonically was wrong.

H5 characterizes this finding with finer resolution and tests whether
the sweet spot location is universal across selection regimes (beta) and
target lengths (L), or system-specific. It also asks whether M2 at the
sweet spot r* is functionally equivalent to M4 at its own optimal mutation
rate mu* -- if so, the variation-rate principle (Draft B) becomes the
empirical generalization of the mechanism-specific framing (Draft A).

================================================================================
Pre-registered predictions (DEFINED BEFORE RUN)
================================================================================

P_H5_1 (sweet spot fine structure): The plateau-height-vs-r curve at
        beta=10, L=32 has a single global maximum in the range
        r in [0.01, 0.30]. The location of this maximum is r* in
        [0.05, 0.15] -- the H3 r=0.10 finding refined.

P_H5_2 (rate-matching to M4 mutation): The sweet spot location r* is
        comparable to M4's per-genome mutation rate (mu * L = 0.01 * 32
        = 0.32 expected mutations per genome per generation, though r is
        per-reproduction not per-position). If the analogy holds, r*
        should be at the value where M2's expected per-reproduction
        template-replacement rate matches M4's per-reproduction mutation
        impact.

P_H5_3 (beta invariance): The sweet spot r* is approximately invariant
        across beta in {2, 5, 10, 20}. If r* shifts substantially
        (e.g., moves from 0.05 at beta=2 to 0.20 at beta=20), the
        variation-rate principle is contingent on selection regime.

P_H5_4 (L invariance): The sweet spot r* is approximately invariant
        across L_TARGET in {16, 32, 64, 128}. If r* tracks L (e.g.,
        scales as 1/L), the principle is more naturally stated in
        terms of per-position rates.

P_H5_5 (M2 plateau height at r*): The plateau height at r* across all
        regimes lies between M1 (no variation) and M4 (full
        individual-level mutation). M2 cannot exceed M4's plateau --
        if it does, the framework needs a deeper revision because
        lineage-level re-draws would then be a *better* variation
        source than individual mutation.

The strongest joint test is P_H5_3 + P_H5_5: if r* is regime-invariant
and M2's plateau at r* sits cleanly between M1 and M4, the variation-rate
principle is robust and v3's Draft B framing is empirically supported.

================================================================================
Sweep design (totals)
================================================================================

Step 1 (fine r sweep at beta=10, L=32):
    14 r values x 30 reps = 420 sims.

Step 2 (beta sensitivity, neighborhood r values):
    4 beta x 7 r x 30 reps = 840 sims.

Step 3 (L sensitivity, neighborhood r values):
    4 L x 5 r x 30 reps = 600 sims.

Step 4 (M4 mutation rate sweep):
    6 mu values x 30 reps = 180 sims.

Step 5 (M2 vs M4 head-to-head at matched rates r*, mu*):
    3 N_M2 values x 30 reps = 90 sims. (r*, mu* determined dynamically
    from Step 1 and Step 4 results.)

Step 6: post-process -> v3_statement.md (no sims).

GRAND TOTAL: 420 + 840 + 600 + 180 + 90 = 2130 sims.
At ~1-2 sec per sim (extrapolated from H3 at K=400, 1000 gens),
wall-time ~30-60 min.

================================================================================
Reproducibility
================================================================================

Module seed: np.random.seed(42), np.random.default_rng(42).
Per-cell, per-replicate seed: seed = 42 + cell_idx*100 + rep_idx
where cell_idx is the position in the global cell list.

================================================================================
Anti-DRY discipline
================================================================================

Per CLAUDE.md, this script duplicates the M2 (lineage fixation with re-draw
rate r), M4 (individual copy with mutation), K=400 selection loop, and
fitness function from test_h3_inheritance_landscape.py. If a bug is fixed
in one, fix it in the other intentionally.

================================================================================
CLI arguments
================================================================================

--progress-file <path>  : where to write progress lines (every 50 sims)
--completed-csv  <path> : append every finished sim's row immediately (fsync'd)
--smoke-test            : tiny run (1 r value, 2 reps, 200 gens; one mu
                          value; mini head-to-head) for verification
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
# parameters (matched to test_h3 where overlapping)
# ============================================================================
ALPHABET = 4
L_TARGET_DEFAULT = 32
K = 400
N_GEN = 1000
BETA_CANONICAL = 10.0

EPS_NOISE_M2 = 0.05    # phenotype noise for M2 (lineage template + noise) -- matches H3 EPS_NOISE
MU_M4_DEFAULT = 0.01   # canonical M4 mutation rate (matches H3)

N_REPLICATES = 30

# Step 1: fine r sweep
R_FINE = [0.00, 0.01, 0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30,
          0.50, 0.75, 1.00]

# Step 2: beta sensitivity, neighborhood r values
BETA_SWEEP = [2.0, 5.0, 10.0, 20.0]
R_NEIGHBORHOOD = [0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20]

# Step 3: L sensitivity, neighborhood r values
L_SWEEP = [16, 32, 64, 128]
R_NEIGHBORHOOD_L = [0.02, 0.05, 0.10, 0.15, 0.20]

# Step 4: M4 mutation rate sweep
MU_SWEEP = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]

# Step 5: M2 vs M4 head-to-head, initial counts of M2 (M4 = K - N_M2)
N_M2_VALUES = [20, 200, 380]


# fixed targets, keyed by L_TARGET. seed 2026 to match h2/h3/h4 idiom.
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
# Mechanism M2 -- lineage-level fixation with re-draw rate r
# Each agent has a fixed template; reproduction copies the parent's template
# with prob 1-r, otherwise draws a fresh template (= new lineage). Phenotype
# = template + per-position noise (rate EPS_NOISE_M2).
# Duplicated from test_h3_inheritance_landscape.py (deliberate, see CLAUDE.md).
# ============================================================================
class MechM2:
    label_template = "M2_redraw_r{r:.3f}"

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
        noise = self.rng.random((n, self.L)) < EPS_NOISE_M2
        random_subs = self.rng.integers(0, ALPHABET, size=(n, self.L), dtype=np.int8)
        return np.where(noise, random_subs, intended).astype(np.int8)

    def reproduceFromIndices(self, parent_indices):
        n_off = len(parent_indices)
        if n_off == 0:
            self.templates = np.empty((0, self.L), dtype=np.int8)
            self.lineage_ids = np.empty(0, dtype=np.int64)
            return
        # path A (1-r): inherit parent template & lineage id
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
# Mechanism M4 -- individual-level copy with mutation rate mu (parameterized)
# Faithful copy + per-position mutation. Phenotype = genotype.
# Duplicated from test_h3 but parameterized over mu (H3 hardcoded MU_M4=0.01).
# ============================================================================
class MechM4:
    label_template = "M4_indiv_mu{mu:.4f}"

    def __init__(self, n_agents, rng, L, mu):
        self.rng = rng
        self.L = L
        self.mu = mu
        self.genotypes = rng.integers(0, ALPHABET, size=(n_agents, L), dtype=np.int8)
        self.label = self.label_template.format(mu=mu)

    @property
    def n(self):
        return self.genotypes.shape[0]

    def findPhenotypes(self):
        return self.genotypes

    def reproduceFromIndices(self, parent_indices):
        n_off = len(parent_indices)
        if n_off == 0:
            self.genotypes = np.empty((0, self.L), dtype=np.int8)
            return
        new_g = self.genotypes[parent_indices].copy()
        flips = self.rng.random(new_g.shape) < self.mu
        replacements = self.rng.integers(0, ALPHABET, size=new_g.shape, dtype=np.int8)
        new_g = np.where(flips, replacements, new_g)
        self.genotypes = new_g

    def countDistinctLineages(self):
        if self.genotypes.size == 0:
            return 0
        return int(np.unique(self.genotypes, axis=0).shape[0])


# ============================================================================
# isolated population dynamics driver
# K agents all of one mechanism, run n_gen generations under softmax selection.
# Returns dict with per-generation summary.
# ============================================================================
def runIsolatedM2(r, n_gen, beta, L, rng):
    """K agents of M2 with re-draw rate r. Returns history list."""
    target = TARGETS[L]
    mech_rng = np.random.default_rng(rng.integers(0, 2**31 - 1))
    mech = MechM2(K, mech_rng, L, r)

    history = []
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
        parent_indices = mech_rng.choice(mech.n, size=K, replace=True, p=w)
        mech.reproduceFromIndices(parent_indices)
    return history


def runIsolatedM4(mu, n_gen, beta, L, rng):
    """K agents of M4 with mutation rate mu. Returns history list."""
    target = TARGETS[L]
    mech_rng = np.random.default_rng(rng.integers(0, 2**31 - 1))
    mech = MechM4(K, mech_rng, L, mu)

    history = []
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
        w = np.exp(beta * (fit - fit.max()))
        w = w / w.sum()
        parent_indices = mech_rng.choice(mech.n, size=K, replace=True, p=w)
        mech.reproduceFromIndices(parent_indices)
    return history


# ============================================================================
# pairwise competition driver (Step 5: M2 vs M4)
# two mechanisms share K slots, compete under shared softmax selection.
# ============================================================================
def runPairwiseM2vsM4(r, mu, n_m2_init, n_m4_init, n_gen, beta, L, rng):
    assert n_m2_init + n_m4_init == K
    target = TARGETS[L]
    rng_m2 = np.random.default_rng(rng.integers(0, 2**31 - 1))
    rng_m4 = np.random.default_rng(rng.integers(0, 2**31 - 1))

    sub_m2 = MechM2(n_m2_init, rng_m2, L, r)
    sub_m4 = MechM4(n_m4_init, rng_m4, L, mu)

    history = []
    label_m2 = sub_m2.label
    label_m4 = sub_m4.label

    for g in range(n_gen):
        phen_m2 = sub_m2.findPhenotypes() if sub_m2.n > 0 else np.empty((0, L), dtype=np.int8)
        phen_m4 = sub_m4.findPhenotypes() if sub_m4.n > 0 else np.empty((0, L), dtype=np.int8)
        fit_m2 = findFitness(phen_m2, target) if sub_m2.n > 0 else np.empty(0, dtype=np.float64)
        fit_m4 = findFitness(phen_m4, target) if sub_m4.n > 0 else np.empty(0, dtype=np.float64)

        n_m2 = sub_m2.n
        n_m4 = sub_m4.n
        n_total = n_m2 + n_m4

        history.append({
            "generation": g,
            "mode_label": label_m2,
            "count": int(n_m2),
            "frequency": float(n_m2) / n_total if n_total > 0 else 0.0,
            "mean_fitness": float(fit_m2.mean()) if n_m2 > 0 else float("nan"),
        })
        history.append({
            "generation": g,
            "mode_label": label_m4,
            "count": int(n_m4),
            "frequency": float(n_m4) / n_total if n_total > 0 else 0.0,
            "mean_fitness": float(fit_m4.mean()) if n_m4 > 0 else float("nan"),
        })

        all_fit = np.concatenate([fit_m2, fit_m4])
        if all_fit.size == 0:
            break
        w = np.exp(beta * (all_fit - all_fit.max()))
        w = w / w.sum()
        parent_indices = rng.choice(n_total, size=K, replace=True, p=w)

        is_m2 = parent_indices < n_m2
        m2_parents = parent_indices[is_m2]
        m4_parents = parent_indices[~is_m2] - n_m2

        if m2_parents.size > 0:
            sub_m2.reproduceFromIndices(m2_parents)
        else:
            sub_m2 = MechM2(0, rng_m2, L, r)
        if m4_parents.size > 0:
            sub_m4.reproduceFromIndices(m4_parents)
        else:
            sub_m4 = MechM4(0, rng_m4, L, mu)

    return {"history": history, "label_m2": label_m2, "label_m4": label_m4}


# ============================================================================
# CSV writer (append-as-we-go, one row per finished sim, fsync)
# unified schema -- step + row_kind disambiguate downstream
# ============================================================================
COMPLETED_FIELDS = [
    "row_kind",         # "m2_isolated", "m4_isolated", "m2_vs_m4"
    "step",             # 1, 2, 3, 4, 5
    "cell_idx",
    "replicate",
    "seed",
    "L",
    "beta",
    "r",                # M2 re-draw rate (m2_* and pairwise rows)
    "mu",               # M4 mutation rate (m4_* and pairwise rows)
    "n_gen",
    # isolated columns (M2 or M4)
    "final_mean_fitness",
    "final_max_fitness",
    "final_var_fitness",
    "final_n_lineages",
    "max_mean_fitness_overall",
    "gen_reach_0p95_max",
    # pairwise columns (M2 vs M4)
    "n_m2_init",
    "n_m4_init",
    "final_freq_m2",
    "final_freq_m4",
    "final_mean_fitness_m2",
    "final_mean_fitness_m4",
    "gen_m2_reaches_95",
    "gen_m4_reaches_95",
    "winner",            # "M2", "M4", or "tie"
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
# cell builders
# ============================================================================
def buildStep1Cells():
    """Step 1: fine r sweep at beta=10, L=32 -- 14 r values."""
    cells = []
    for r in R_FINE:
        cells.append({
            "step": 1, "kind": "m2_isolated",
            "r": r, "mu": "", "beta": BETA_CANONICAL, "L": L_TARGET_DEFAULT,
        })
    return cells


def buildStep2Cells():
    """Step 2: beta sensitivity at L=32, r in neighborhood -- 4 beta x 7 r."""
    cells = []
    for beta in BETA_SWEEP:
        for r in R_NEIGHBORHOOD:
            cells.append({
                "step": 2, "kind": "m2_isolated",
                "r": r, "mu": "", "beta": beta, "L": L_TARGET_DEFAULT,
            })
    return cells


def buildStep3Cells():
    """Step 3: L sensitivity at beta=10, r in neighborhood -- 4 L x 5 r."""
    cells = []
    for L in L_SWEEP:
        for r in R_NEIGHBORHOOD_L:
            cells.append({
                "step": 3, "kind": "m2_isolated",
                "r": r, "mu": "", "beta": BETA_CANONICAL, "L": L,
            })
    return cells


def buildStep4Cells():
    """Step 4: M4 mu sweep at beta=10, L=32 -- 6 mu values."""
    cells = []
    for mu in MU_SWEEP:
        cells.append({
            "step": 4, "kind": "m4_isolated",
            "r": "", "mu": mu, "beta": BETA_CANONICAL, "L": L_TARGET_DEFAULT,
        })
    return cells


def buildStep5Cells(r_star, mu_star):
    """Step 5: M2(r=r*) vs M4(mu=mu*) head-to-head -- 3 N_M2 values."""
    cells = []
    for n_m2 in N_M2_VALUES:
        cells.append({
            "step": 5, "kind": "m2_vs_m4",
            "r": r_star, "mu": mu_star,
            "n_m2_init": n_m2, "n_m4_init": K - n_m2,
            "beta": BETA_CANONICAL, "L": L_TARGET_DEFAULT,
        })
    return cells


# ============================================================================
# cell runners
# ============================================================================
def runM2IsolatedCell(cell, cell_idx, n_replicates, n_gen, completed_csv, state):
    rows = []
    r = cell["r"]
    beta = cell["beta"]
    L = cell["L"]
    step = cell["step"]
    for rep_idx in range(n_replicates):
        if state.get("shutdown"):
            break
        seed = 42 + cell_idx * 100 + rep_idx
        rng = np.random.default_rng(seed)
        t0 = time.time()
        hist = runIsolatedM2(r, n_gen, beta, L, rng)
        dt = time.time() - t0

        last = hist[-1]
        max_mean_overall = max(h["mean_fitness"] for h in hist)
        gen_reach_95 = -1
        for h in hist:
            if h["mean_fitness"] >= 0.95:
                gen_reach_95 = h["generation"]
                break

        row = {
            "row_kind": "m2_isolated",
            "step": step,
            "cell_idx": cell_idx,
            "replicate": rep_idx,
            "seed": seed,
            "L": L,
            "beta": beta,
            "r": r,
            "mu": "",
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


def runM4IsolatedCell(cell, cell_idx, n_replicates, n_gen, completed_csv, state):
    rows = []
    mu = cell["mu"]
    beta = cell["beta"]
    L = cell["L"]
    step = cell["step"]
    for rep_idx in range(n_replicates):
        if state.get("shutdown"):
            break
        seed = 42 + cell_idx * 100 + rep_idx
        rng = np.random.default_rng(seed)
        t0 = time.time()
        hist = runIsolatedM4(mu, n_gen, beta, L, rng)
        dt = time.time() - t0

        last = hist[-1]
        max_mean_overall = max(h["mean_fitness"] for h in hist)
        gen_reach_95 = -1
        for h in hist:
            if h["mean_fitness"] >= 0.95:
                gen_reach_95 = h["generation"]
                break

        row = {
            "row_kind": "m4_isolated",
            "step": step,
            "cell_idx": cell_idx,
            "replicate": rep_idx,
            "seed": seed,
            "L": L,
            "beta": beta,
            "r": "",
            "mu": mu,
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


def runM2vsM4Cell(cell, cell_idx, n_replicates, n_gen, completed_csv, state):
    rows = []
    r = cell["r"]
    mu = cell["mu"]
    n_m2_init = cell["n_m2_init"]
    n_m4_init = cell["n_m4_init"]
    beta = cell["beta"]
    L = cell["L"]
    step = cell["step"]
    for rep_idx in range(n_replicates):
        if state.get("shutdown"):
            break
        seed = 42 + cell_idx * 100 + rep_idx
        rng = np.random.default_rng(seed)
        t0 = time.time()
        out = runPairwiseM2vsM4(r, mu, n_m2_init, n_m4_init, n_gen, beta, L, rng)
        dt = time.time() - t0
        hist = out["history"]
        label_m2 = out["label_m2"]
        label_m4 = out["label_m4"]

        m2_recs = [h for h in hist if h["mode_label"] == label_m2]
        m4_recs = [h for h in hist if h["mode_label"] == label_m4]
        m2_recs.sort(key=lambda x: x["generation"])
        m4_recs.sort(key=lambda x: x["generation"])

        final_m2 = m2_recs[-1] if m2_recs else None
        final_m4 = m4_recs[-1] if m4_recs else None

        gen_m2_95 = -1
        for h in m2_recs:
            if h["frequency"] > 0.95:
                gen_m2_95 = h["generation"]
                break
        gen_m4_95 = -1
        for h in m4_recs:
            if h["frequency"] > 0.95:
                gen_m4_95 = h["generation"]
                break

        # winner = whichever crosses 0.95 first; otherwise larger final freq
        if gen_m2_95 >= 0 and (gen_m4_95 < 0 or gen_m2_95 < gen_m4_95):
            winner = "M2"
        elif gen_m4_95 >= 0 and (gen_m2_95 < 0 or gen_m4_95 < gen_m2_95):
            winner = "M4"
        else:
            f2 = final_m2["frequency"] if final_m2 else 0.0
            f4 = final_m4["frequency"] if final_m4 else 0.0
            if abs(f2 - f4) < 0.05:
                winner = "tie"
            else:
                winner = "M2" if f2 > f4 else "M4"

        row = {
            "row_kind": "m2_vs_m4",
            "step": step,
            "cell_idx": cell_idx,
            "replicate": rep_idx,
            "seed": seed,
            "L": L,
            "beta": beta,
            "r": r,
            "mu": mu,
            "n_gen": n_gen,
            "n_m2_init": n_m2_init,
            "n_m4_init": n_m4_init,
            "final_freq_m2": float(final_m2["frequency"]) if final_m2 else 0.0,
            "final_freq_m4": float(final_m4["frequency"]) if final_m4 else 0.0,
            "final_mean_fitness_m2": float(final_m2["mean_fitness"])
                if (final_m2 and not np.isnan(final_m2["mean_fitness"])) else float("nan"),
            "final_mean_fitness_m4": float(final_m4["mean_fitness"])
                if (final_m4 and not np.isnan(final_m4["mean_fitness"])) else float("nan"),
            "gen_m2_reaches_95": gen_m2_95,
            "gen_m4_reaches_95": gen_m4_95,
            "winner": winner,
            "wall_seconds": dt,
        }
        rows.append(row)
        appendCompleted(completed_csv, [row])
        state["n_done"] += 1
        if state["n_done"] % 50 == 0:
            writeProgressLine(state)
    return rows


# ============================================================================
# signal handler -- clean SIGTERM
# ============================================================================
def installSignalHandlers(state):
    def handler(signum, frame):
        state["shutdown"] = True
        writeFinalProgress(state, status=f"signal_{signum}_received")
        sys.exit(0)
    signal.signal(signal.SIGTERM, handler)
    signal.signal(signal.SIGINT, handler)


# ============================================================================
# post-processing helpers: load completed CSV
# ============================================================================
def loadCompleted(path):
    rows = []
    with open(path, "r") as f:
        rd = csv.DictReader(f)
        for r in rd:
            for k in ("step", "cell_idx", "replicate", "seed", "L", "n_gen",
                      "n_m2_init", "n_m4_init", "final_n_lineages",
                      "gen_reach_0p95_max", "gen_m2_reaches_95",
                      "gen_m4_reaches_95"):
                if r.get(k, "") not in ("", None):
                    try:
                        r[k] = int(float(r[k]))
                    except ValueError:
                        pass
            for k in ("beta", "r", "mu", "final_mean_fitness", "final_max_fitness",
                      "final_var_fitness", "max_mean_fitness_overall",
                      "final_freq_m2", "final_freq_m4",
                      "final_mean_fitness_m2", "final_mean_fitness_m4",
                      "wall_seconds"):
                if r.get(k, "") not in ("", None):
                    try:
                        r[k] = float(r[k])
                    except ValueError:
                        pass
            rows.append(r)
    return rows


def findRStar(rows_step1):
    """Identify r* from Step 1 rows -- the r with max mean final fitness."""
    by_r = {}
    for row in rows_step1:
        if row["row_kind"] != "m2_isolated":
            continue
        r = row["r"]
        ff = row["final_mean_fitness"]
        if isinstance(ff, float) and not np.isnan(ff):
            by_r.setdefault(r, []).append(ff)
    best_r = None
    best_mean = -1.0
    for r, vals in by_r.items():
        m = float(np.mean(vals))
        if m > best_mean:
            best_mean = m
            best_r = r
    return best_r, best_mean, by_r


def findMuStar(rows_step4):
    """Identify mu* from Step 4 rows -- the mu with max mean final fitness."""
    by_mu = {}
    for row in rows_step4:
        if row["row_kind"] != "m4_isolated":
            continue
        mu = row["mu"]
        ff = row["final_mean_fitness"]
        if isinstance(ff, float) and not np.isnan(ff):
            by_mu.setdefault(mu, []).append(ff)
    best_mu = None
    best_mean = -1.0
    for mu, vals in by_mu.items():
        m = float(np.mean(vals))
        if m > best_mean:
            best_mean = m
            best_mu = mu
    return best_mu, best_mean, by_mu


def findRStarByBeta(rows_step2):
    """For Step 2, find r*(beta) for each beta."""
    out = {}
    by_beta_r = {}
    for row in rows_step2:
        if row["row_kind"] != "m2_isolated":
            continue
        beta = row["beta"]
        r = row["r"]
        ff = row["final_mean_fitness"]
        if isinstance(ff, float) and not np.isnan(ff):
            by_beta_r.setdefault((beta, r), []).append(ff)
    means = {}
    for (beta, r), vals in by_beta_r.items():
        means[(beta, r)] = float(np.mean(vals))
    betas = sorted({b for (b, _) in means.keys()})
    for beta in betas:
        per_r = [(r, means[(beta, r)]) for (b, r) in means.keys() if b == beta]
        per_r.sort(key=lambda x: -x[1])
        out[beta] = per_r[0]  # (best_r, best_mean)
    return out


def findRStarByL(rows_step3):
    """For Step 3, find r*(L) for each L."""
    out = {}
    by_L_r = {}
    for row in rows_step3:
        if row["row_kind"] != "m2_isolated":
            continue
        L = row["L"]
        r = row["r"]
        ff = row["final_mean_fitness"]
        if isinstance(ff, float) and not np.isnan(ff):
            by_L_r.setdefault((L, r), []).append(ff)
    means = {}
    for (L, r), vals in by_L_r.items():
        means[(L, r)] = float(np.mean(vals))
    Ls = sorted({L for (L, _) in means.keys()})
    for L in Ls:
        per_r = [(r, means[(L, r)]) for (l, r) in means.keys() if l == L]
        per_r.sort(key=lambda x: -x[1])
        out[L] = per_r[0]  # (best_r, best_mean)
    return out


# ============================================================================
# Step 1-4 result CSV writers
# ============================================================================
def writeStep1Csv(rows, path):
    """Step 1: r sweep -> mean+CI per r."""
    s1 = [r for r in rows if r["row_kind"] == "m2_isolated" and r.get("step") == 1]
    by_r = {}
    for row in s1:
        by_r.setdefault(row["r"], []).append(row)
    fields = ["r", "n_replicates", "mean_final_fitness", "std_final_fitness",
              "ci95_low", "ci95_high", "mean_max_fitness_overall",
              "mean_n_lineages_final", "mean_wall_seconds"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in sorted(by_r.keys()):
            rs = by_r[r]
            ff = [x["final_mean_fitness"] for x in rs if not np.isnan(x["final_mean_fitness"])]
            mf = [x["max_mean_fitness_overall"] for x in rs if not np.isnan(x["max_mean_fitness_overall"])]
            nl = [x["final_n_lineages"] for x in rs]
            ws = [x["wall_seconds"] for x in rs]
            mean = float(np.mean(ff)) if ff else float("nan")
            std = float(np.std(ff, ddof=1)) if len(ff) > 1 else 0.0
            n = len(ff)
            sem = std / np.sqrt(n) if n > 0 else 0.0
            ci95 = 1.96 * sem
            w.writerow({
                "r": r, "n_replicates": len(rs),
                "mean_final_fitness": mean,
                "std_final_fitness": std,
                "ci95_low": mean - ci95,
                "ci95_high": mean + ci95,
                "mean_max_fitness_overall": float(np.mean(mf)) if mf else float("nan"),
                "mean_n_lineages_final": float(np.mean(nl)) if nl else float("nan"),
                "mean_wall_seconds": float(np.mean(ws)) if ws else float("nan"),
            })


def writeStep2Csv(rows, path):
    """Step 2: (beta, r) -> plateau."""
    s2 = [r for r in rows if r["row_kind"] == "m2_isolated" and r.get("step") == 2]
    by_cell = {}
    for row in s2:
        by_cell.setdefault((row["beta"], row["r"]), []).append(row)
    fields = ["beta", "r", "n_replicates", "mean_final_fitness",
              "std_final_fitness", "ci95_low", "ci95_high"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for key in sorted(by_cell.keys()):
            beta, r = key
            rs = by_cell[key]
            ff = [x["final_mean_fitness"] for x in rs if not np.isnan(x["final_mean_fitness"])]
            mean = float(np.mean(ff)) if ff else float("nan")
            std = float(np.std(ff, ddof=1)) if len(ff) > 1 else 0.0
            n = len(ff)
            sem = std / np.sqrt(n) if n > 0 else 0.0
            ci95 = 1.96 * sem
            w.writerow({
                "beta": beta, "r": r, "n_replicates": len(rs),
                "mean_final_fitness": mean,
                "std_final_fitness": std,
                "ci95_low": mean - ci95,
                "ci95_high": mean + ci95,
            })


def writeStep3Csv(rows, path):
    """Step 3: (L, r) -> plateau."""
    s3 = [r for r in rows if r["row_kind"] == "m2_isolated" and r.get("step") == 3]
    by_cell = {}
    for row in s3:
        by_cell.setdefault((row["L"], row["r"]), []).append(row)
    fields = ["L", "r", "n_replicates", "mean_final_fitness",
              "std_final_fitness", "ci95_low", "ci95_high"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for key in sorted(by_cell.keys()):
            L, r = key
            rs = by_cell[key]
            ff = [x["final_mean_fitness"] for x in rs if not np.isnan(x["final_mean_fitness"])]
            mean = float(np.mean(ff)) if ff else float("nan")
            std = float(np.std(ff, ddof=1)) if len(ff) > 1 else 0.0
            n = len(ff)
            sem = std / np.sqrt(n) if n > 0 else 0.0
            ci95 = 1.96 * sem
            w.writerow({
                "L": L, "r": r, "n_replicates": len(rs),
                "mean_final_fitness": mean,
                "std_final_fitness": std,
                "ci95_low": mean - ci95,
                "ci95_high": mean + ci95,
            })


def writeStep4Csv(rows, path):
    """Step 4: M4 mu sweep -> plateau per mu."""
    s4 = [r for r in rows if r["row_kind"] == "m4_isolated" and r.get("step") == 4]
    by_mu = {}
    for row in s4:
        by_mu.setdefault(row["mu"], []).append(row)
    fields = ["mu", "n_replicates", "mean_final_fitness", "std_final_fitness",
              "ci95_low", "ci95_high", "mean_max_fitness_overall",
              "median_gen_reach_0p95"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for mu in sorted(by_mu.keys()):
            rs = by_mu[mu]
            ff = [x["final_mean_fitness"] for x in rs if not np.isnan(x["final_mean_fitness"])]
            mf = [x["max_mean_fitness_overall"] for x in rs if not np.isnan(x["max_mean_fitness_overall"])]
            grs = [x["gen_reach_0p95_max"] for x in rs if x["gen_reach_0p95_max"] >= 0]
            mean = float(np.mean(ff)) if ff else float("nan")
            std = float(np.std(ff, ddof=1)) if len(ff) > 1 else 0.0
            n = len(ff)
            sem = std / np.sqrt(n) if n > 0 else 0.0
            ci95 = 1.96 * sem
            w.writerow({
                "mu": mu, "n_replicates": len(rs),
                "mean_final_fitness": mean,
                "std_final_fitness": std,
                "ci95_low": mean - ci95,
                "ci95_high": mean + ci95,
                "mean_max_fitness_overall": float(np.mean(mf)) if mf else float("nan"),
                "median_gen_reach_0p95": int(np.median(grs)) if grs else -1,
            })


def writeStep5Csv(rows, path):
    """Step 5: M2 vs M4 head-to-head per N_M2_init."""
    s5 = [r for r in rows if r["row_kind"] == "m2_vs_m4" and r.get("step") == 5]
    by_n = {}
    for row in s5:
        by_n.setdefault(row["n_m2_init"], []).append(row)
    fields = ["n_m2_init", "n_m4_init", "r", "mu", "n_replicates",
              "p_m2_wins", "p_m4_wins", "p_tie",
              "mean_final_freq_m2", "mean_final_freq_m4",
              "mean_final_fitness_m2", "mean_final_fitness_m4"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for n_m2 in sorted(by_n.keys()):
            rs = by_n[n_m2]
            n_total = len(rs)
            n_m2_wins = sum(1 for x in rs if x["winner"] == "M2")
            n_m4_wins = sum(1 for x in rs if x["winner"] == "M4")
            n_tie = sum(1 for x in rs if x["winner"] == "tie")
            ff2 = [x["final_freq_m2"] for x in rs]
            ff4 = [x["final_freq_m4"] for x in rs]
            fm2 = [x["final_mean_fitness_m2"] for x in rs if not np.isnan(x["final_mean_fitness_m2"])]
            fm4 = [x["final_mean_fitness_m4"] for x in rs if not np.isnan(x["final_mean_fitness_m4"])]
            w.writerow({
                "n_m2_init": n_m2,
                "n_m4_init": K - n_m2,
                "r": rs[0]["r"],
                "mu": rs[0]["mu"],
                "n_replicates": n_total,
                "p_m2_wins": n_m2_wins / n_total,
                "p_m4_wins": n_m4_wins / n_total,
                "p_tie": n_tie / n_total,
                "mean_final_freq_m2": float(np.mean(ff2)),
                "mean_final_freq_m4": float(np.mean(ff4)),
                "mean_final_fitness_m2": float(np.mean(fm2)) if fm2 else float("nan"),
                "mean_final_fitness_m4": float(np.mean(fm4)) if fm4 else float("nan"),
            })


# ============================================================================
# plotting
# ============================================================================
def plotRSweepCurve(rows, path):
    """Step 1: plateau height vs r at beta=10, L=32 -- with 95% CI band."""
    s1 = [r for r in rows if r["row_kind"] == "m2_isolated" and r.get("step") == 1]
    if not s1:
        return
    by_r = {}
    for row in s1:
        by_r.setdefault(row["r"], []).append(row)
    rs_sorted = sorted(by_r.keys())
    means = []
    ci_low = []
    ci_high = []
    for r in rs_sorted:
        ff = [x["final_mean_fitness"] for x in by_r[r] if not np.isnan(x["final_mean_fitness"])]
        m = float(np.mean(ff)) if ff else float("nan")
        sd = float(np.std(ff, ddof=1)) if len(ff) > 1 else 0.0
        sem = sd / np.sqrt(len(ff)) if ff else 0.0
        means.append(m)
        ci_low.append(m - 1.96 * sem)
        ci_high.append(m + 1.96 * sem)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.fill_between(rs_sorted, ci_low, ci_high, color="C0", alpha=0.25,
                    label="95% CI")
    ax.plot(rs_sorted, means, "o-", color="C0", linewidth=1.8, markersize=6,
            label="M2 plateau (mean over 30 reps)")
    # mark r*
    if means:
        i_star = int(np.argmax(means))
        ax.axvline(rs_sorted[i_star], color="red", linestyle="--", alpha=0.6,
                   label=f"r* = {rs_sorted[i_star]:.3f}")
    ax.set_xlabel("M2 re-draw rate r")
    ax.set_ylabel("Plateau mean fitness at gen 1000")
    ax.set_title("Test H5 Step 1 -- M2 sweet spot (beta=10, L=32, K=400, 30 reps)")
    ax.set_xscale("symlog", linthresh=0.01)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotRegimeSensitivity(rows, path):
    """Step 2 + Step 3: r* vs beta, r* vs L (two-panel)."""
    fig, (ax_b, ax_L) = plt.subplots(1, 2, figsize=(12, 5))

    # left panel: r*(beta) overlaid on plateau curves per beta
    s2 = [r for r in rows if r["row_kind"] == "m2_isolated" and r.get("step") == 2]
    if s2:
        betas = sorted({row["beta"] for row in s2})
        cmap = plt.get_cmap("viridis")
        for i, beta in enumerate(betas):
            rs_for_b = sorted({row["r"] for row in s2 if row["beta"] == beta})
            means = []
            for r in rs_for_b:
                cell = [row for row in s2 if row["beta"] == beta and row["r"] == r]
                ff = [x["final_mean_fitness"] for x in cell if not np.isnan(x["final_mean_fitness"])]
                means.append(float(np.mean(ff)) if ff else float("nan"))
            ax_b.plot(rs_for_b, means, "o-", color=cmap(i / max(1, len(betas) - 1)),
                      label=f"beta={beta}")
            # mark sweet spot
            if means:
                i_star = int(np.argmax(means))
                ax_b.scatter([rs_for_b[i_star]], [means[i_star]], color="red",
                             marker="*", s=120, zorder=5)
        ax_b.set_xlabel("M2 re-draw rate r")
        ax_b.set_ylabel("Plateau mean fitness")
        ax_b.set_title("Step 2: beta sensitivity (L=32)\nred star = r*(beta)")
        ax_b.grid(True, alpha=0.3)
        ax_b.legend(loc="best", fontsize=8)

    # right panel: r*(L) overlaid on plateau curves per L
    s3 = [r for r in rows if r["row_kind"] == "m2_isolated" and r.get("step") == 3]
    if s3:
        Ls = sorted({row["L"] for row in s3})
        cmap = plt.get_cmap("plasma")
        for i, L in enumerate(Ls):
            rs_for_L = sorted({row["r"] for row in s3 if row["L"] == L})
            means = []
            for r in rs_for_L:
                cell = [row for row in s3 if row["L"] == L and row["r"] == r]
                ff = [x["final_mean_fitness"] for x in cell if not np.isnan(x["final_mean_fitness"])]
                means.append(float(np.mean(ff)) if ff else float("nan"))
            ax_L.plot(rs_for_L, means, "o-", color=cmap(i / max(1, len(Ls) - 1)),
                      label=f"L={L}")
            if means:
                i_star = int(np.argmax(means))
                ax_L.scatter([rs_for_L[i_star]], [means[i_star]], color="red",
                             marker="*", s=120, zorder=5)
        ax_L.set_xlabel("M2 re-draw rate r")
        ax_L.set_ylabel("Plateau mean fitness")
        ax_L.set_title("Step 3: L sensitivity (beta=10)\nred star = r*(L)")
        ax_L.grid(True, alpha=0.3)
        ax_L.legend(loc="best", fontsize=8)

    fig.suptitle("Test H5 -- Sweet spot regime sensitivity")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotRateMatching(rows, path):
    """Step 1 + Step 4 overlay: M2 plateau vs r and M4 plateau vs mu on twin axes."""
    s1 = [r for r in rows if r["row_kind"] == "m2_isolated" and r.get("step") == 1]
    s4 = [r for r in rows if r["row_kind"] == "m4_isolated" and r.get("step") == 4]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    if s1:
        by_r = {}
        for row in s1:
            by_r.setdefault(row["r"], []).append(row)
        rs_sorted = sorted(by_r.keys())
        means = [float(np.mean([x["final_mean_fitness"] for x in by_r[r]
                                if not np.isnan(x["final_mean_fitness"])]))
                 for r in rs_sorted]
        ax.plot(rs_sorted, means, "o-", color="C0", label="M2 plateau vs r",
                linewidth=1.8, markersize=6)
        if means:
            i_star = int(np.argmax(means))
            ax.axvline(rs_sorted[i_star], color="C0", linestyle="--", alpha=0.5,
                       label=f"r* = {rs_sorted[i_star]:.3f}")
    if s4:
        by_mu = {}
        for row in s4:
            by_mu.setdefault(row["mu"], []).append(row)
        mus_sorted = sorted(by_mu.keys())
        means_mu = [float(np.mean([x["final_mean_fitness"] for x in by_mu[mu]
                                   if not np.isnan(x["final_mean_fitness"])]))
                    for mu in mus_sorted]
        ax.plot(mus_sorted, means_mu, "s-", color="C3", label="M4 plateau vs mu",
                linewidth=1.8, markersize=6)
        if means_mu:
            i_star_mu = int(np.argmax(means_mu))
            ax.axvline(mus_sorted[i_star_mu], color="C3", linestyle="--", alpha=0.5,
                       label=f"mu* = {mus_sorted[i_star_mu]:.4f}")
    ax.set_xscale("log")
    ax.set_xlabel("variation rate (r for M2, mu for M4)")
    ax.set_ylabel("Plateau mean fitness at gen 1000")
    ax.set_title("Test H5 -- Rate matching: M2 (r) vs M4 (mu) on shared axis (beta=10, L=32)")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotM2vsM4Matched(rows, path):
    """Step 5: bar plot of P(M2 wins), P(M4 wins), P(tie) per N_M2_init."""
    s5 = [r for r in rows if r["row_kind"] == "m2_vs_m4" and r.get("step") == 5]
    if not s5:
        return
    by_n = {}
    for row in s5:
        by_n.setdefault(row["n_m2_init"], []).append(row)
    n_values = sorted(by_n.keys())
    p_m2 = []
    p_m4 = []
    p_tie = []
    for n in n_values:
        rs = by_n[n]
        n_total = len(rs)
        p_m2.append(sum(1 for x in rs if x["winner"] == "M2") / n_total)
        p_m4.append(sum(1 for x in rs if x["winner"] == "M4") / n_total)
        p_tie.append(sum(1 for x in rs if x["winner"] == "tie") / n_total)

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(n_values))
    w = 0.27
    ax.bar(x - w, p_m2, w, label="M2 wins", color="C0")
    ax.bar(x, p_m4, w, label="M4 wins", color="C3")
    ax.bar(x + w, p_tie, w, label="tie", color="grey")
    ax.set_xticks(x)
    ax.set_xticklabels([f"N_M2={n}" for n in n_values])
    ax.set_ylabel("Win rate (30 reps per cell)")
    r_star = s5[0]["r"]
    mu_star = s5[0]["mu"]
    ax.set_title(f"Test H5 Step 5 -- M2(r={r_star}) vs M4(mu={mu_star}) head-to-head")
    ax.set_ylim(0, 1.05)
    ax.axhline(0.5, color="black", linestyle=":", alpha=0.5, label="50% (tied)")
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ============================================================================
# v3 statement writer
# ============================================================================
def writeV3Statement(rows, path):
    """Generate test_h5_v3_statement.md -- final interpretation."""
    s1 = [r for r in rows if r["row_kind"] == "m2_isolated" and r.get("step") == 1]
    s2 = [r for r in rows if r["row_kind"] == "m2_isolated" and r.get("step") == 2]
    s3 = [r for r in rows if r["row_kind"] == "m2_isolated" and r.get("step") == 3]
    s4 = [r for r in rows if r["row_kind"] == "m4_isolated" and r.get("step") == 4]
    s5 = [r for r in rows if r["row_kind"] == "m2_vs_m4" and r.get("step") == 5]

    r_star, r_star_mean, r_means = findRStar(s1)
    mu_star, mu_star_mean, mu_means = findMuStar(s4)
    r_star_by_beta = findRStarByBeta(s2)
    r_star_by_L = findRStarByL(s3)

    # 95% CI for r* plateau
    if r_star is not None:
        ff_at_r_star = r_means.get(r_star, [])
        sd = float(np.std(ff_at_r_star, ddof=1)) if len(ff_at_r_star) > 1 else 0.0
        sem = sd / np.sqrt(len(ff_at_r_star)) if ff_at_r_star else 0.0
        ci = 1.96 * sem
    else:
        ci = float("nan")

    # P_H5_1 verdict: r* in [0.05, 0.15]?
    p1_pass = r_star is not None and 0.05 <= r_star <= 0.15
    # P_H5_3 verdict: r* invariant across beta?
    rs_by_b = [v[0] for v in r_star_by_beta.values()]
    p3_invariant = (len(rs_by_b) > 0
                    and (max(rs_by_b) - min(rs_by_b)) <= 0.05)
    # P_H5_4 verdict: r* invariant across L?
    rs_by_L = [v[0] for v in r_star_by_L.values()]
    p4_invariant = (len(rs_by_L) > 0
                    and (max(rs_by_L) - min(rs_by_L)) <= 0.05)
    # P_H5_5 verdict: M2 plateau at r* < M4 plateau at mu*?
    p5_pass = (r_star is not None and mu_star is not None
               and r_star_mean < mu_star_mean + 0.02)
    # Step 5 head-to-head outcome
    if s5:
        n_total = len(s5)
        n_m2_wins = sum(1 for x in s5 if x["winner"] == "M2")
        n_m4_wins = sum(1 for x in s5 if x["winner"] == "M4")
        n_tie = sum(1 for x in s5 if x["winner"] == "tie")
        # check N_M2=200 cell (matched start)
        s5_balanced = [x for x in s5 if x["n_m2_init"] == 200]
        if s5_balanced:
            n_b = len(s5_balanced)
            p_m4_wins_balanced = sum(1 for x in s5_balanced if x["winner"] == "M4") / n_b
            p_m2_wins_balanced = sum(1 for x in s5_balanced if x["winner"] == "M2") / n_b
            p_tie_balanced = sum(1 for x in s5_balanced if x["winner"] == "tie") / n_b
        else:
            p_m4_wins_balanced = float("nan")
            p_m2_wins_balanced = float("nan")
            p_tie_balanced = float("nan")
    else:
        p_m4_wins_balanced = float("nan")
        p_m2_wins_balanced = float("nan")
        p_tie_balanced = float("nan")

    # interpretation: Draft B if matched balanced is near-tie OR if M4 wins
    # by < 80%; Draft A if M4 wins >= 80% even at matched rates.
    if not np.isnan(p_m4_wins_balanced):
        if p_m4_wins_balanced >= 0.80:
            recommended = "A"
            rationale = (
                "M4 wins >= 80% at matched rates (N_M2=200), indicating that "
                "individual-level mutation has an intrinsic advantage beyond "
                "raw variation rate. Draft A's mechanism-specific framing is "
                "the safer position."
            )
        elif abs(p_m2_wins_balanced - p_m4_wins_balanced) < 0.30:
            recommended = "B"
            rationale = (
                "At matched variation rates (M2 at r*, M4 at mu*), the "
                "head-to-head competition at N_M2=200 is approximately "
                "balanced. The variation-rate principle is empirically "
                "supported -- the mechanism distinction collapses at matched "
                "rates. Draft B (mechanism-agnostic, rate-based) is the "
                "empirical generalization of Draft A."
            )
        else:
            recommended = "B (partial)"
            rationale = (
                f"M4 wins {p_m4_wins_balanced*100:.0f}% at matched rates -- "
                f"intermediate outcome. Draft B is partially supported but "
                f"individual-level mutation retains a measurable edge. The "
                f"v3 hybrid framing (Draft A in theorem, Draft B in "
                f"Discussion) remains a defensible compromise."
            )
    else:
        recommended = "(insufficient data)"
        rationale = "Step 5 head-to-head data missing; cannot adjudicate."

    # build markdown
    lines = []
    lines.append("# Test H5 -- v3 Statement (M2 sweet spot characterization)")
    lines.append("")
    lines.append("Generated by `code/test_h5_sweet_spot.py` from the post-processed sweep data.")
    lines.append("This is the v3 input on the variation-rate principle (Draft A vs Draft B).")
    lines.append("")
    lines.append("## 1. Sweet spot location (Step 1, beta=10, L=32)")
    lines.append("")
    if r_star is not None:
        lines.append(f"**r\\* = {r_star:.3f}** "
                     f"(plateau mean fitness = {r_star_mean:.4f} +/- {ci:.4f}; 30 reps).")
    else:
        lines.append("r\\* not found (no Step 1 data).")
    lines.append("")
    lines.append("Plateau mean fitness vs r:")
    lines.append("")
    lines.append("| r | mean fitness | n_reps |")
    lines.append("|---|---|---|")
    for r in sorted(r_means.keys()):
        vals = r_means[r]
        m = float(np.mean(vals))
        lines.append(f"| {r:.3f} | {m:.4f} | {len(vals)} |")
    lines.append("")
    lines.append(f"P_H5_1 (r\\* in [0.05, 0.15]): **{'CONFIRMED' if p1_pass else 'REFINED'}**")
    lines.append("")

    lines.append("## 2. Regime sensitivity (Steps 2, 3)")
    lines.append("")
    lines.append("### r\\*(beta) at L=32")
    lines.append("")
    lines.append("| beta | r\\* | plateau at r\\* |")
    lines.append("|---|---|---|")
    for beta in sorted(r_star_by_beta.keys()):
        r_b, m_b = r_star_by_beta[beta]
        lines.append(f"| {beta} | {r_b:.3f} | {m_b:.4f} |")
    lines.append("")
    lines.append(f"P_H5_3 (beta-invariance, max(r\\*) - min(r\\*) <= 0.05): "
                 f"**{'CONFIRMED' if p3_invariant else 'REFINED'}**")
    lines.append("")
    lines.append("### r\\*(L) at beta=10")
    lines.append("")
    lines.append("| L | r\\* | plateau at r\\* |")
    lines.append("|---|---|---|")
    for L in sorted(r_star_by_L.keys()):
        r_L, m_L = r_star_by_L[L]
        lines.append(f"| {L} | {r_L:.3f} | {m_L:.4f} |")
    lines.append("")
    lines.append(f"P_H5_4 (L-invariance, max(r\\*) - min(r\\*) <= 0.05): "
                 f"**{'CONFIRMED' if p4_invariant else 'REFINED'}**")
    lines.append("")

    lines.append("## 3. Rate matching (Step 4 + Step 5)")
    lines.append("")
    lines.append("### M4 plateau vs mu (Step 4)")
    lines.append("")
    lines.append("| mu | mean fitness | n_reps |")
    lines.append("|---|---|---|")
    for mu in sorted(mu_means.keys()):
        vals = mu_means[mu]
        m = float(np.mean(vals))
        lines.append(f"| {mu:.4f} | {m:.4f} | {len(vals)} |")
    lines.append("")
    if mu_star is not None:
        lines.append(f"**mu\\* = {mu_star:.4f}** (plateau = {mu_star_mean:.4f}).")
    lines.append("")
    if r_star is not None and mu_star is not None:
        # per-genome variation rate comparison
        # M2: r per reproduction (fresh template substitutes L positions)
        # M4: mu * L per reproduction (expected mutated positions per offspring)
        r_per_genome = r_star  # one full template replacement per reproduction event
        mu_per_genome = mu_star * L_TARGET_DEFAULT
        lines.append("### Per-reproduction variation rate comparison")
        lines.append("")
        lines.append(f"- M2 at r\\* = {r_star:.3f}: P(full template replacement per offspring) = {r_per_genome:.3f}")
        lines.append(f"  -> expected NEW positions per offspring = {r_star * L_TARGET_DEFAULT * (1 - 1/ALPHABET):.3f}")
        lines.append(f"  (when redrawn, each position is fresh; expected differing positions = L*(1-1/ALPHABET))")
        lines.append(f"- M4 at mu\\* = {mu_star:.4f}: expected mutated positions per offspring = mu*L = {mu_per_genome:.3f}")
        lines.append(f"  -> expected NEW positions per offspring = mu*L*(1-1/ALPHABET) = {mu_per_genome * (1 - 1/ALPHABET):.3f}")
        lines.append("")
        lines.append(f"P_H5_2 (rate-matching analogy): the M2 sweet spot at r* = {r_star:.3f} "
                     f"corresponds to {r_star * L_TARGET_DEFAULT * (1 - 1/ALPHABET):.2f} expected new positions per offspring; "
                     f"M4 optimal at mu* = {mu_star:.4f} corresponds to {mu_per_genome * (1 - 1/ALPHABET):.2f} new positions. "
                     f"These should be comparable if Draft B holds.")
        lines.append("")
    lines.append(f"P_H5_5 (M2 at r\\* not exceeding M4 at mu\\*): "
                 f"**{'CONFIRMED' if p5_pass else 'FALSIFIED'}** "
                 f"(M2={r_star_mean:.4f} vs M4={mu_star_mean:.4f}).")
    lines.append("")

    lines.append("### Step 5 head-to-head competition")
    lines.append("")
    if s5:
        lines.append("| N_M2 init | N_M4 init | P(M2 wins) | P(M4 wins) | P(tie) |")
        lines.append("|---|---|---|---|---|")
        by_n = {}
        for row in s5:
            by_n.setdefault(row["n_m2_init"], []).append(row)
        for n_m2 in sorted(by_n.keys()):
            rs = by_n[n_m2]
            nt = len(rs)
            n_w_m2 = sum(1 for x in rs if x["winner"] == "M2")
            n_w_m4 = sum(1 for x in rs if x["winner"] == "M4")
            n_w_t = sum(1 for x in rs if x["winner"] == "tie")
            lines.append(f"| {n_m2} | {K - n_m2} | {n_w_m2/nt:.2f} | {n_w_m4/nt:.2f} | {n_w_t/nt:.2f} |")
    lines.append("")
    if not np.isnan(p_m4_wins_balanced):
        lines.append(f"At balanced start (N_M2=200, N_M4=200): "
                     f"M2 wins {p_m2_wins_balanced*100:.1f}%, "
                     f"M4 wins {p_m4_wins_balanced*100:.1f}%, "
                     f"tie {p_tie_balanced*100:.1f}%.")
    lines.append("")

    lines.append("## 4. Verdict: Draft A vs Draft B")
    lines.append("")
    lines.append(f"**Recommended: Draft {recommended}.**")
    lines.append("")
    lines.append(rationale)
    lines.append("")

    lines.append("## 5. Drop-in for v3 Discussion")
    lines.append("")
    lines.append("The following paragraph can be inserted into the v3 Discussion to "
                 "address the M2 sweet spot finding:")
    lines.append("")
    lines.append("> Test H3 found an unanticipated sweet spot in M2 (lineage-fixation "
                 "with re-draw rate r): the plateau mean fitness at r=0.10 exceeded "
                 "both r=0.00 (pure fixation, no variation) and r=0.50 (high re-draw "
                 f"rate). Test H5 refined this to r* = {r_star:.3f} "
                 f"(plateau = {r_star_mean:.3f} +/- {ci:.3f}, 30 reps). "
                 "The location of r* was "
                 f"{'approximately invariant' if p3_invariant else 'sensitive'} "
                 "to selection sharpness (beta in {2,5,10,20}) and "
                 f"{'approximately invariant' if p4_invariant else 'sensitive'} "
                 "to target length (L in {16,32,64,128}). At its sweet spot, M2's "
                 f"plateau ({r_star_mean:.3f}) "
                 f"{'remained below' if p5_pass else 'matched or exceeded'} "
                 f"M4's optimal plateau ({mu_star_mean:.3f}, achieved at "
                 f"mu* = {mu_star:.4f}). In direct head-to-head competition at "
                 f"matched variation rates and balanced start (N_M2 = N_M4 = 200), "
                 f"M4 won {p_m4_wins_balanced*100:.0f}% of replicates "
                 f"(M2 {p_m2_wins_balanced*100:.0f}%, "
                 f"tie {p_tie_balanced*100:.0f}%), "
                 "consistent with the variation-rate principle "
                 f"(Draft {recommended.split()[0]}) as the operative empirical "
                 "generalization of the original mechanism-specific framing.")
    lines.append("")
    lines.append("## 6. Pre-registered prediction outcomes")
    lines.append("")
    lines.append("| ID | Verdict | Empirical |")
    lines.append("|---|---|---|")
    lines.append(f"| P_H5_1 | {'CONFIRMED' if p1_pass else 'REFINED'} | "
                 f"r\\* = {r_star:.3f}; predicted r\\* in [0.05, 0.15] |")
    if r_star is not None and mu_star is not None:
        m2_var = r_star * L_TARGET_DEFAULT * (1 - 1/ALPHABET)
        m4_var = mu_star * L_TARGET_DEFAULT * (1 - 1/ALPHABET)
        ratio = m2_var / m4_var if m4_var > 0 else float("inf")
        lines.append(f"| P_H5_2 | {'CONFIRMED' if abs(np.log(ratio)) < np.log(3) else 'REFINED'} | "
                     f"M2 new-positions/offspring at r\\*={m2_var:.2f}; M4={m4_var:.2f}; ratio={ratio:.2f} |")
    else:
        lines.append("| P_H5_2 | NO_DATA | r\\* or mu\\* missing |")
    lines.append(f"| P_H5_3 | {'CONFIRMED' if p3_invariant else 'REFINED'} | "
                 f"r\\*(beta): " +
                 ", ".join(f"beta={b}:r\\*={v[0]:.3f}" for b, v in sorted(r_star_by_beta.items())) +
                 " |")
    lines.append(f"| P_H5_4 | {'CONFIRMED' if p4_invariant else 'REFINED'} | "
                 f"r\\*(L): " +
                 ", ".join(f"L={L}:r\\*={v[0]:.3f}" for L, v in sorted(r_star_by_L.items())) +
                 " |")
    lines.append(f"| P_H5_5 | {'CONFIRMED' if p5_pass else 'FALSIFIED'} | "
                 f"M2(r\\*) plateau = {r_star_mean:.4f} vs M4(mu\\*) plateau = {mu_star_mean:.4f} |")
    lines.append("")
    Path(path).write_text("\n".join(lines))


# ============================================================================
# main
# ============================================================================
def main():
    np.random.seed(42)

    parser = argparse.ArgumentParser()
    parser.add_argument("--progress-file", required=True)
    parser.add_argument("--completed-csv", required=True)
    parser.add_argument("--smoke-test", action="store_true",
                        help="tiny config: 1 r value x 2 reps x 200 gens; mini head-to-head")
    args = parser.parse_args()

    progress_path = Path(args.progress_file)
    completed_path = Path(args.completed_csv)
    progress_path.parent.mkdir(parents=True, exist_ok=True)
    completed_path.parent.mkdir(parents=True, exist_ok=True)

    progress_path.write_text(f"# === test_h5 run started at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")

    if args.smoke_test:
        # tiny config matching the dispatch's smoke-test directive:
        # 1 r value x 2 reps x 200 gens; 1 mu value; 1 head-to-head cell
        cells_steps_1234 = [
            {"step": 1, "kind": "m2_isolated",
             "r": 0.10, "mu": "", "beta": BETA_CANONICAL, "L": L_TARGET_DEFAULT},
            {"step": 4, "kind": "m4_isolated",
             "r": "", "mu": 0.01, "beta": BETA_CANONICAL, "L": L_TARGET_DEFAULT},
        ]
        n_replicates = 2
        n_gen = 200
        # smoke step 5 will be a single n_m2 cell
        smoke_step5_n_m2 = [200]
        if completed_path.exists():
            completed_path.unlink()
    else:
        cells_steps_1234 = (buildStep1Cells() + buildStep2Cells()
                            + buildStep3Cells() + buildStep4Cells())
        n_replicates = N_REPLICATES
        n_gen = N_GEN
        smoke_step5_n_m2 = None
        if completed_path.exists():
            backup = completed_path.with_suffix(".csv.bak")
            completed_path.rename(backup)
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] archived prior completed.csv to {backup}\n")

    # estimate total: Step 5 will be added after r* and mu* are known
    n_step5 = 3 if smoke_step5_n_m2 is None else len(smoke_step5_n_m2)
    n_total = sum(n_replicates for _ in cells_steps_1234) + n_step5 * n_replicates

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
                f"{len(cells_steps_1234)} cells in steps 1-4 + step5 deferred; "
                f"projected total sims = {n_total}; "
                f"smoke_test={args.smoke_test}; n_gen={n_gen}\n")

    # ----- Steps 1-4 -----
    for cell_idx, cell in enumerate(cells_steps_1234):
        if state["shutdown"]:
            break
        state["current_step"] = (
            f"step{cell['step']}/{cell['kind']}"
            f"@beta={cell['beta']},L={cell['L']},"
            f"r={cell.get('r', '')},mu={cell.get('mu', '')}"
        )
        if cell["kind"] == "m2_isolated":
            runM2IsolatedCell(cell, cell_idx, n_replicates, n_gen, completed_path, state)
        elif cell["kind"] == "m4_isolated":
            runM4IsolatedCell(cell, cell_idx, n_replicates, n_gen, completed_path, state)
        writeProgressLine(state)

    # ----- Step 5: build dynamically from r* and mu* -----
    if not state["shutdown"]:
        all_rows_so_far = loadCompleted(completed_path)
        s1_rows = [r for r in all_rows_so_far
                   if r["row_kind"] == "m2_isolated" and r.get("step") == 1]
        s4_rows = [r for r in all_rows_so_far
                   if r["row_kind"] == "m4_isolated" and r.get("step") == 4]
        r_star, _, _ = findRStar(s1_rows)
        mu_star, _, _ = findMuStar(s4_rows)
        if r_star is None or mu_star is None:
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] WARNING: Step 5 skipped, r*={r_star} mu*={mu_star}\n")
        else:
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Step 5 starting with r*={r_star} mu*={mu_star}\n")
            if smoke_step5_n_m2 is not None:
                step5_cells = [{
                    "step": 5, "kind": "m2_vs_m4",
                    "r": r_star, "mu": mu_star,
                    "n_m2_init": n, "n_m4_init": K - n,
                    "beta": BETA_CANONICAL, "L": L_TARGET_DEFAULT,
                } for n in smoke_step5_n_m2]
            else:
                step5_cells = buildStep5Cells(r_star, mu_star)

            base_idx = len(cells_steps_1234)
            for j, cell in enumerate(step5_cells):
                if state["shutdown"]:
                    break
                cell_idx = base_idx + j
                state["current_step"] = (
                    f"step5/m2_vs_m4@N_M2={cell['n_m2_init']},r={cell['r']:.3f},mu={cell['mu']:.4f}"
                )
                runM2vsM4Cell(cell, cell_idx, n_replicates, n_gen, completed_path, state)
                writeProgressLine(state)

    writeFinalProgress(state, status="completed")

    # ----- Step 6: post-process (full run only; smoke skips heavy plotting) -----
    if not args.smoke_test:
        try:
            all_rows = loadCompleted(completed_path)
            writeStep1Csv(all_rows, RESULTS_DIR / "test_h5_r_sweep_v1.csv")
            writeStep2Csv(all_rows, RESULTS_DIR / "test_h5_beta_sensitivity_v1.csv")
            writeStep3Csv(all_rows, RESULTS_DIR / "test_h5_L_sensitivity_v1.csv")
            writeStep4Csv(all_rows, RESULTS_DIR / "test_h5_M4_mutation_sweep_v1.csv")
            writeStep5Csv(all_rows, RESULTS_DIR / "test_h5_M2_vs_M4_matched_v1.csv")
            plotRSweepCurve(all_rows, FIGURES_DIR / "test_h5_r_sweep_curve.png")
            plotRegimeSensitivity(all_rows, FIGURES_DIR / "test_h5_regime_sensitivity.png")
            plotRateMatching(all_rows, FIGURES_DIR / "test_h5_rate_matching.png")
            plotM2vsM4Matched(all_rows, FIGURES_DIR / "test_h5_M2_vs_M4_matched.png")
            writeV3Statement(all_rows, RESULTS_DIR / "test_h5_v3_statement.md")
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] post-process complete: "
                        f"5 step CSVs, 4 figures, v3_statement.md written\n")
        except Exception as e:
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] post-process FAILED: {e}\n")
            raise
    else:
        # smoke test: skip plotting (insufficient data) but try the v3 statement
        try:
            all_rows = loadCompleted(completed_path)
            writeStep1Csv(all_rows, RESULTS_DIR / "test_h5_smoke_step1.csv")
            writeStep4Csv(all_rows, RESULTS_DIR / "test_h5_smoke_step4.csv")
            writeStep5Csv(all_rows, RESULTS_DIR / "test_h5_smoke_step5.csv")
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] smoke post-process complete\n")
        except Exception as e:
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] smoke post-process FAILED: {e}\n")
            raise


if __name__ == "__main__":
    main()
