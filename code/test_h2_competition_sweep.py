"""
Test H2 -- Founder-loss boundary + graded copyability sweep.

================================================================================
Purpose
================================================================================

Test H v1 found P_H3 (the strict abundance-independence claim) failed in
1/10 reps of scenario B: Mode 1 went stochastically extinct from 40 founders
(Wright-Fisher founder loss). H2 maps the founder-loss boundary across
initial Mode 1 abundances N_1 and tests whether *graded* copyability (h in
[0,1]) for Mode 6 changes the outcome. The result supersedes the strict
P_H3 statement of H v1 with a quantitatively-bounded claim suitable for
the v3 paper Discussion.

================================================================================
Pre-registered predictions (DEFINED BEFORE RUN)
================================================================================

P_H2_1 (founder-loss threshold): For Mode 1 vs Mode 6 (Implementation A)
        competition with Mode 1 initial counts N_1 in {20,40,60,80,100,120,
        160,200,240,280,320,360,380} out of K=400 and 30 reps per cell, the
        probability of Mode 1 reaching >95% by gen 1000 increases monotonically
        with N_1, transitioning from <50% at N_1 = 20 to >95% at N_1 >= 80.

P_H2_2 (Implementation A robustness): The founder-loss boundary identified in
        P_H2_1 holds for Implementation A (binary copyability, fixed template
        per replicate) -- the canonical case in v1.

P_H2_3 (Implementation B robustness): For Implementation B (per-lineage fixed
        template, no inheritance of mutations), Mode 6 should behave
        identically to Implementation A because neither implements copy. The
        boundary should match A's at all N_1 values, within replicate noise.

P_H2_4 (Implementation C graded copyability): For Implementation C
        (heritability h in {0.0, 0.25, 0.5, 0.75, 1.0}), Mode 6 should
        increasingly resemble Mode 1 as h -> 1. At h = 0.5, Mode 6 should
        compete with Mode 1 but still lose because Mode 1's full heritability
        gives it a fitness-climb advantage. At h = 1.0, Mode 6 IS Mode 1
        (the framework's claim is binary, so this should be a degenerate case).

P_H2_5 (graded crossover): The crossover N_1 (initial Mode 1 count where
        P(Mode 1 wins) crosses 0.5) should shift right as Mode 6's heritability
        h increases. At h = 0, the crossover N_1 is the founder-loss threshold
        from P_H2_1. At h = 1, the crossover is at N_1 = 200 (the 50/50
        starting point) by symmetry.

The strongest test is P_H2_5: if h = 0.5 doesn't shift the crossover
significantly relative to h = 0, then framework's "binary copyability matters"
claim is challenged by the data.

================================================================================
Implementations of Mode 6 (parameterization rule)
================================================================================

Implementation A (binary copyability, replicate-fixed template):
    Mode 6 has a SINGLE fixed template shared by the entire Mode 6 population
    for the entire replicate. Each agent each generation draws phenotype =
    template + i.i.d. noise (rate EPS_NOISE_6). No copy mechanism: offspring
    just get a fresh draw. Identical to test_h v1.

Implementation B (per-lineage fixed template, no mutation inheritance):
    Each lineage gets its own random template at founding; the template is
    fixed for that lineage's entire history and is INHERITED by descendants.
    Mutations to phenotype noise do not propagate. New lineages (at the start)
    get distinct random templates. The framework's prediction (P_H2_3): B
    should match A because neither has copy mechanism.

Implementation C (graded copyability, heritability h):
    Mode 6 reproduction: with probability h, offspring's genotype = parent's
    genotype (faithful copy with mutation MU_MODE1, like Mode 1). With
    probability (1-h), offspring's genotype is replaced by a fresh draw from
    the fixed-template + noise distribution (like Implementation A). At h=0
    this is Implementation A; at h=1 it is Mode 1 with the same mutation rate.

================================================================================
Sweep design (totals)
================================================================================

Step 1 (Implementation A): N_1 in {20,40,60,80,100,120,160,200,240,280,320,
        360,380} = 13 cells x 30 reps = 390 sims
Step 2 (Implementation B): N_1 in {20,40,80,120,160,200} = 6 cells x 30 reps
        = 180 sims
Step 3 (Implementation C): h in {0.0,0.25,0.5,0.75,1.0} x N_1 in
        {20,40,80,120,160,200,240,280,320,360,380} = 5 x 11 = 55 cells
        x 30 reps = 1650 sims
Total: 390 + 180 + 1650 = 2220 simulations
Each sim: K=400, N_GEN=1000, ~5-15 sec on Boron (vectorized over 400 agents).

================================================================================
Reproducibility
================================================================================

Module-level seed: np.random.seed(42), np.random.default_rng(42).
Per-cell, per-replicate seed scheme:
    seed = 42 + cell_idx * 100 + rep_idx
where cell_idx is a global integer assigned by enumeration order across all
(impl, h, N_1) combinations and rep_idx in {0..29}. This gives each (cell, rep)
a unique seed and avoids collisions across the three implementations.

================================================================================
Parameters (matched to test_h v1 / test_d_v2)
================================================================================

ALPHABET    = 4
L_TARGET    = 32
K           = 400
N_GEN       = 1000
BETA        = 10.0
MU_MODE1    = 0.01
N_MODULES_5 = 8 (unused in H2; kept for symmetry)
EPS_NOISE_5 = 0.05 (unused)
EPS_NOISE_6 = 0.05
N_REPLICATES = 30 (per cell)

================================================================================
Anti-DRY discipline
================================================================================

Per CLAUDE.md, this script duplicates the Mode 1 / Mode 6 logic from
test_h_competition.py and test_d_v2_population_dynamics.py rather than
importing. If you fix a bug in one, fix it intentionally in the others.

================================================================================
CLI arguments
================================================================================

--progress-file <path>  : where to write progress lines (every 50 sims)
--completed-csv  <path>  : append every finished sim's row immediately
                            (recovery-friendly)
--smoke-test            : tiny test (Implementation A only, 2 N_1 values
                            x 2 reps x 50 generations) for verification
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
# parameters (must match test_h v1 / test_d_v2)
# ============================================================================
ALPHABET = 4
L_TARGET = 32
K = 400
N_GEN = 1000
BETA = 10.0

MU_MODE1 = 0.01
EPS_NOISE_6 = 0.05

N_REPLICATES = 30

# fixed target shared across the whole H2 run, same idiom as test_h v1.
TARGET = np.random.default_rng(2026).integers(0, ALPHABET, size=L_TARGET, dtype=np.int8)


# ============================================================================
# fitness (vectorized)
# ============================================================================
def findFitness(phenotypes):
    """phenotypes shape (N, L_TARGET) -> fitness (N,) in [0,1]."""
    matches = (phenotypes == TARGET[None, :]).sum(axis=1)
    return matches.astype(np.float64) / L_TARGET


# ============================================================================
# Mode 1: 1D template, faithful copy with point mutation MU_MODE1
# (duplicated from test_h_competition.py; do not import)
# ============================================================================
class Mode1Subpop:
    label = "Mode1"

    def __init__(self, n_agents, rng):
        self.rng = rng
        self.templates = rng.integers(0, ALPHABET, size=(n_agents, L_TARGET), dtype=np.int8)

    @property
    def n(self):
        return self.templates.shape[0]

    def findPhenotypes(self):
        # mode 1 phenotype = template (no noise per test_d_v2 / test_h v1)
        return self.templates

    def reproduceFromIndices(self, parent_indices):
        new_templates = self.templates[parent_indices].copy()
        flips = self.rng.random(new_templates.shape) < MU_MODE1
        replacements = self.rng.integers(0, ALPHABET, size=new_templates.shape, dtype=np.int8)
        new_templates = np.where(flips, replacements, new_templates)
        self.templates = new_templates


# ============================================================================
# Mode 6 -- Implementation A (replicate-fixed template, no copy)
# duplicated from test_h_competition.py
# ============================================================================
class Mode6SubpopA:
    label = "Mode6_implA"

    def __init__(self, n_agents, rng, fixed_template):
        self.rng = rng
        self.fixed_template = fixed_template  # shape (L_TARGET,)
        self._n = n_agents

    @property
    def n(self):
        return self._n

    def findPhenotypes(self):
        n_agents = self._n
        if n_agents == 0:
            return np.empty((0, L_TARGET), dtype=np.int8)
        intended = np.broadcast_to(self.fixed_template[None, :], (n_agents, L_TARGET))
        noise = self.rng.random((n_agents, L_TARGET)) < EPS_NOISE_6
        random_subs = self.rng.integers(0, ALPHABET, size=(n_agents, L_TARGET), dtype=np.int8)
        return np.where(noise, random_subs, intended).astype(np.int8)

    def reproduceFromIndices(self, parent_indices):
        # no genotype, no copying; only the count carries forward.
        self._n = len(parent_indices)


# ============================================================================
# Mode 6 -- Implementation B (per-lineage fixed template, INHERITED to descendants)
# Each agent has a fixed template; offspring inherits parent's template.
# Phenotype draw = template + i.i.d. noise. No mutation in template.
# Initial founders get DISTINCT random templates (one per founder).
# ============================================================================
class Mode6SubpopB:
    label = "Mode6_implB"

    def __init__(self, n_agents, rng):
        self.rng = rng
        # per-agent template; founders get distinct random templates
        self.templates = rng.integers(0, ALPHABET, size=(n_agents, L_TARGET), dtype=np.int8)

    @property
    def n(self):
        return self.templates.shape[0]

    def findPhenotypes(self):
        n_agents = self.templates.shape[0]
        if n_agents == 0:
            return np.empty((0, L_TARGET), dtype=np.int8)
        # phenotype = own template + i.i.d. noise
        intended = self.templates
        noise = self.rng.random((n_agents, L_TARGET)) < EPS_NOISE_6
        random_subs = self.rng.integers(0, ALPHABET, size=(n_agents, L_TARGET), dtype=np.int8)
        return np.where(noise, random_subs, intended).astype(np.int8)

    def reproduceFromIndices(self, parent_indices):
        # offspring inherits parent's template exactly; no template mutation.
        # this models lineage-fixed templates with no copy mechanism --
        # adaptation cannot accumulate within a lineage.
        self.templates = self.templates[parent_indices].copy()


# ============================================================================
# Mode 6 -- Implementation C (graded copyability with heritability h)
# With probability h, offspring inherits parent genotype (with Mode 1
# mutation). With probability (1-h), offspring genotype is overwritten by
# a fresh draw from the fixed-template + noise distribution.
# Genotype is the per-agent template; phenotype = genotype directly (h=1
# becomes Mode 1 exactly).
# At h = 0: behaves like Implementation A (template-fixed noise).
# At h = 1: behaves like Mode 1 with the same mutation rate.
# ============================================================================
class Mode6SubpopC:
    label_template = "Mode6_implC_h{h:.2f}"

    def __init__(self, n_agents, rng, fixed_template, h):
        self.rng = rng
        self.fixed_template = fixed_template  # shape (L_TARGET,)
        self.h = h
        # initialize all agents with a fresh draw from the fixed-template
        # + noise distribution (same as the "non-inherit" path).
        self.genotypes = self._drawFresh(n_agents)
        self.label = self.label_template.format(h=h)

    @property
    def n(self):
        return self.genotypes.shape[0]

    def _drawFresh(self, n_agents):
        # fresh sample from fixed template + i.i.d. noise (rate EPS_NOISE_6).
        # this is the (1 - h) path; matches Implementation A's noise model.
        if n_agents == 0:
            return np.empty((0, L_TARGET), dtype=np.int8)
        intended = np.broadcast_to(self.fixed_template[None, :], (n_agents, L_TARGET))
        noise = self.rng.random((n_agents, L_TARGET)) < EPS_NOISE_6
        random_subs = self.rng.integers(0, ALPHABET, size=(n_agents, L_TARGET), dtype=np.int8)
        return np.where(noise, random_subs, intended).astype(np.int8)

    def findPhenotypes(self):
        # phenotype = genotype directly (so at h=1, this is Mode 1)
        return self.genotypes

    def reproduceFromIndices(self, parent_indices):
        n_off = len(parent_indices)
        if n_off == 0:
            self.genotypes = np.empty((0, L_TARGET), dtype=np.int8)
            return
        # path 1 (with prob h): inherit parent genotype + Mode 1 mutation
        copied = self.genotypes[parent_indices].copy()
        flips = self.rng.random(copied.shape) < MU_MODE1
        replacements = self.rng.integers(0, ALPHABET, size=copied.shape, dtype=np.int8)
        copied = np.where(flips, replacements, copied)
        # path 2 (with prob 1-h): fresh draw from fixed-template + noise
        fresh = self._drawFresh(n_off)
        # per-agent draw to decide which path each offspring takes
        use_copy = self.rng.random(n_off) < self.h
        self.genotypes = np.where(use_copy[:, None], copied, fresh)


# ============================================================================
# competition driver: two subpops sharing K slots
# duplicated from test_h_competition.py
# ============================================================================
def runCompetition(make_subpop_a, make_subpop_b, n_a_init, n_b_init, n_gen, rng):
    """Coupled two-mode competition for n_gen generations.

    make_subpop_a / make_subpop_b: callables (n_agents, rng) -> Subpop
    n_a_init + n_b_init must equal K.
    Returns: list of per-generation per-mode dicts.
    """
    assert n_a_init + n_b_init == K, f"initial counts must sum to K={K}"

    rng_a = np.random.default_rng(rng.integers(0, 2**31 - 1))
    rng_b = np.random.default_rng(rng.integers(0, 2**31 - 1))

    sub_a = make_subpop_a(n_a_init, rng_a)
    sub_b = make_subpop_b(n_b_init, rng_b)

    history = []

    for g in range(n_gen):
        # 1. each subpop produces phenotypes for this generation
        phen_a = sub_a.findPhenotypes() if sub_a.n > 0 else np.empty((0, L_TARGET), dtype=np.int8)
        phen_b = sub_b.findPhenotypes() if sub_b.n > 0 else np.empty((0, L_TARGET), dtype=np.int8)

        fit_a = findFitness(phen_a) if sub_a.n > 0 else np.empty(0, dtype=np.float64)
        fit_b = findFitness(phen_b) if sub_b.n > 0 else np.empty(0, dtype=np.float64)

        n_a = sub_a.n
        n_b = sub_b.n
        n_total = n_a + n_b
        history.append({
            "generation": g,
            "mode_label": sub_a.label,
            "count": int(n_a),
            "frequency": float(n_a) / n_total if n_total > 0 else 0.0,
            "mean_fitness": float(fit_a.mean()) if n_a > 0 else float("nan"),
            "max_fitness": float(fit_a.max()) if n_a > 0 else float("nan"),
            "total_population": int(n_total),
        })
        history.append({
            "generation": g,
            "mode_label": sub_b.label,
            "count": int(n_b),
            "frequency": float(n_b) / n_total if n_total > 0 else 0.0,
            "mean_fitness": float(fit_b.mean()) if n_b > 0 else float("nan"),
            "max_fitness": float(fit_b.max()) if n_b > 0 else float("nan"),
            "total_population": int(n_total),
        })

        # 2. shared selection: stack fitnesses, exp(BETA * fit), draw K parents
        all_fit = np.concatenate([fit_a, fit_b])
        if all_fit.size == 0:
            break
        w = np.exp(BETA * (all_fit - all_fit.max()))
        w = w / w.sum()
        parent_indices = rng.choice(n_total, size=K, replace=True, p=w)

        # 3. partition parents back to subpops
        is_a = parent_indices < n_a
        a_parents = parent_indices[is_a]
        b_parents = parent_indices[~is_a] - n_a

        if a_parents.size > 0:
            sub_a.reproduceFromIndices(a_parents)
        else:
            sub_a = make_subpop_a(0, rng_a)
        if b_parents.size > 0:
            sub_b.reproduceFromIndices(b_parents)
        else:
            sub_b = make_subpop_b(0, rng_b)

    return history


# ============================================================================
# per-sim summary helpers
# ============================================================================
def summarizeSim(history, mode_a_label, mode_b_label):
    """Reduce history to one row per mode with key stats."""
    rows = {mode_a_label: [], mode_b_label: []}
    for r in history:
        if r["mode_label"] in rows:
            rows[r["mode_label"]].append(r)

    summary = {}
    for mode in (mode_a_label, mode_b_label):
        rs = sorted(rows[mode], key=lambda x: x["generation"])
        if not rs:
            summary[mode] = {
                "final_count": 0,
                "final_frequency": 0.0,
                "final_mean_fitness": float("nan"),
                "max_mean_fitness": float("nan"),
                "gen_extinct_lt_5pct": -1,
                "gen_dom_gt_95pct": -1,
                "gen_extinct_zero_within_50": -1,
            }
            continue
        gens = [r["generation"] for r in rs]
        freqs = [r["frequency"] for r in rs]
        counts = [r["count"] for r in rs]
        fits = [r["mean_fitness"] for r in rs]
        ext_gen = -1
        for g, f in zip(gens, freqs):
            if f < 0.05:
                ext_gen = g
                break
        dom_gen = -1
        for g, f in zip(gens, freqs):
            if f > 0.95:
                dom_gen = g
                break
        ext_zero_50 = -1
        for g, c in zip(gens, counts):
            if g > 50:
                break
            if c == 0:
                ext_zero_50 = g
                break
        finite_fits = [x for x in fits if not (x is None or np.isnan(x))]
        summary[mode] = {
            "final_count": int(rs[-1]["count"]),
            "final_frequency": float(rs[-1]["frequency"]),
            "final_mean_fitness": float(rs[-1]["mean_fitness"]) if not np.isnan(rs[-1]["mean_fitness"]) else float("nan"),
            "max_mean_fitness": max(finite_fits) if finite_fits else float("nan"),
            "gen_extinct_lt_5pct": ext_gen,
            "gen_dom_gt_95pct": dom_gen,
            "gen_extinct_zero_within_50": ext_zero_50,
        }
    return summary


# ============================================================================
# CSV output: append-as-we-go (one row per (cell, rep, mode) pair)
# ============================================================================
COMPLETED_FIELDS = [
    "implementation", "h", "n1_init", "n6_init", "cell_idx", "replicate", "seed",
    "mode_label",
    "final_count", "final_frequency", "final_mean_fitness", "max_mean_fitness",
    "gen_extinct_lt_5pct", "gen_dom_gt_95pct", "gen_extinct_zero_within_50",
    "wall_seconds",
]


def appendCompleted(path, rows):
    """Append rows to completed.csv; write header if new file."""
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
# implementation factory makers
# ============================================================================
def makeMode1(n_agents, rng):
    return Mode1Subpop(n_agents, rng)


def makeMode6A_factory(fixed_template):
    def maker(n_agents, rng):
        return Mode6SubpopA(n_agents, rng, fixed_template)
    return maker


def makeMode6B(n_agents, rng):
    return Mode6SubpopB(n_agents, rng)


def makeMode6C_factory(fixed_template, h):
    def maker(n_agents, rng):
        return Mode6SubpopC(n_agents, rng, fixed_template, h)
    return maker


# ============================================================================
# cell runner: one (impl, h, N_1) cell over n_replicates
# ============================================================================
def runCell(impl, h, n1_init, cell_idx, n_replicates, n_gen, completed_csv,
            global_state):
    """Runs one (impl, h, N_1) cell of n_replicates sims.

    Appends finished sim rows to completed_csv as it goes.
    Updates global_state['n_done'] for the progress writer.
    Returns the list of summary rows for this cell.
    """
    cell_rows = []
    n6_init = K - n1_init
    label6_for_summary = None  # filled per-rep based on impl

    for rep_idx in range(n_replicates):
        if global_state.get("shutdown"):
            break
        seed = 42 + cell_idx * 100 + rep_idx
        rng = np.random.default_rng(seed)

        # set up Mode 6 maker per implementation
        if impl == "A":
            tmpl_rng = np.random.default_rng(seed * 7919 + 1)
            fixed_template = tmpl_rng.integers(0, ALPHABET, size=L_TARGET, dtype=np.int8)
            make_b = makeMode6A_factory(fixed_template)
            label6 = "Mode6_implA"
        elif impl == "B":
            make_b = makeMode6B
            label6 = "Mode6_implB"
        elif impl == "C":
            tmpl_rng = np.random.default_rng(seed * 7919 + 1)
            fixed_template = tmpl_rng.integers(0, ALPHABET, size=L_TARGET, dtype=np.int8)
            make_b = makeMode6C_factory(fixed_template, h)
            label6 = f"Mode6_implC_h{h:.2f}"
        else:
            raise ValueError(f"unknown impl: {impl}")

        t0 = time.time()
        hist = runCompetition(makeMode1, make_b, n1_init, n6_init, n_gen, rng)
        dt = time.time() - t0

        summary = summarizeSim(hist, "Mode1", label6)

        for mode in ("Mode1", label6):
            row = {
                "implementation": impl,
                "h": h if impl == "C" else "",
                "n1_init": n1_init,
                "n6_init": n6_init,
                "cell_idx": cell_idx,
                "replicate": rep_idx,
                "seed": seed,
                "mode_label": mode,
                "wall_seconds": dt,
                **summary[mode],
            }
            cell_rows.append(row)

        # append the two rows for this sim immediately (recovery-friendly)
        appendCompleted(completed_csv, cell_rows[-2:])
        global_state["n_done"] += 1

        # progress write every 50 sims
        if global_state["n_done"] % 50 == 0:
            writeProgressLine(global_state)

    return cell_rows


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
# sweep configuration
# ============================================================================
STEP1_N1_VALUES = [20, 40, 60, 80, 100, 120, 160, 200, 240, 280, 320, 360, 380]
STEP2_N1_VALUES = [20, 40, 80, 120, 160, 200]
STEP3_H_VALUES = [0.0, 0.25, 0.5, 0.75, 1.0]
STEP3_N1_VALUES = [20, 40, 80, 120, 160, 200, 240, 280, 320, 360, 380]


def buildCellList():
    """Return list of (impl, h, n1) tuples; cell_idx is the position in the list."""
    cells = []
    for n1 in STEP1_N1_VALUES:
        cells.append(("A", 0.0, n1))
    for n1 in STEP2_N1_VALUES:
        cells.append(("B", 0.0, n1))
    for h in STEP3_H_VALUES:
        for n1 in STEP3_N1_VALUES:
            cells.append(("C", h, n1))
    return cells


# ============================================================================
# Step output writers (run after sweep completes from completed.csv)
# ============================================================================
def loadCompletedRows(path):
    rows = []
    with open(path, "r") as f:
        rd = csv.DictReader(f)
        for r in rd:
            # cast numeric fields
            r["n1_init"] = int(r["n1_init"])
            r["n6_init"] = int(r["n6_init"])
            r["cell_idx"] = int(r["cell_idx"])
            r["replicate"] = int(r["replicate"])
            r["seed"] = int(r["seed"])
            for k in ("final_count", "gen_extinct_lt_5pct", "gen_dom_gt_95pct",
                      "gen_extinct_zero_within_50"):
                r[k] = int(float(r[k])) if r[k] != "" else -1
            for k in ("final_frequency", "final_mean_fitness", "max_mean_fitness",
                      "wall_seconds"):
                r[k] = float(r[k]) if r[k] not in ("", "nan") else float("nan")
            r["h"] = float(r["h"]) if r["h"] not in ("",) else 0.0
            rows.append(r)
    return rows


def writeStepCsvs(all_rows):
    """Write per-step result CSVs from the completed-rows DB."""
    # implementation A
    a_rows = [r for r in all_rows if r["implementation"] == "A"]
    writePerCellCsv(a_rows, RESULTS_DIR / "test_h2_founder_boundary_v1.csv",
                    include_h=False)

    # implementation B
    b_rows = [r for r in all_rows if r["implementation"] == "B"]
    writePerCellCsv(b_rows, RESULTS_DIR / "test_h2_implementation_B_v1.csv",
                    include_h=False)

    # implementation C
    c_rows = [r for r in all_rows if r["implementation"] == "C"]
    writePerCellCsv(c_rows, RESULTS_DIR / "test_h2_implementation_C_v1.csv",
                    include_h=True)


def writePerCellCsv(rows, path, include_h):
    """Aggregate per-cell stats and write a CSV."""
    cells = {}
    for r in rows:
        key = (r["implementation"], r["h"], r["n1_init"])
        cells.setdefault(key, []).append(r)

    fields = ["implementation", "n1_init", "n6_init", "n_replicates",
              "p_mode1_wins_gt_95",  # P(Mode 1 reaches >95%)
              "p_mode1_founder_extinct_within_50",
              "median_crossover_gen",  # median gen Mode 1 reached >95% (across reps that won)
              "mean_final_fitness_mode1_winning",
              "mean_final_fitness_mode1_losing",
              "mean_final_freq_mode1",
              "sd_final_freq_mode1",
              "ci_lo_p_wins", "ci_hi_p_wins"]
    if include_h:
        fields.insert(1, "h")

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for key, cell_rows in sorted(cells.items()):
            impl, h, n1 = key
            mode1_rows = [r for r in cell_rows if r["mode_label"] == "Mode1"]
            n_reps = len(mode1_rows)
            if n_reps == 0:
                continue
            # P(Mode 1 wins >95%)
            n_wins = sum(1 for r in mode1_rows if r["final_frequency"] > 0.95)
            p_wins = n_wins / n_reps
            # 95% CI via Wilson score (works for small samples)
            ci_lo, ci_hi = wilsonCi(n_wins, n_reps)
            # P(Mode 1 founder extinct within 50 gens)
            n_founder_ext = sum(1 for r in mode1_rows
                                if r["gen_extinct_zero_within_50"] >= 0)
            p_founder_ext = n_founder_ext / n_reps
            # median crossover gen across reps that won
            won_gens = [r["gen_dom_gt_95pct"] for r in mode1_rows
                        if r["gen_dom_gt_95pct"] >= 0]
            median_cross = int(np.median(won_gens)) if won_gens else -1
            # final fitness conditional
            won_fits = [r["final_mean_fitness"] for r in mode1_rows
                        if r["final_frequency"] > 0.95
                        and not np.isnan(r["final_mean_fitness"])]
            lost_fits = [r["final_mean_fitness"] for r in mode1_rows
                         if r["final_frequency"] <= 0.95
                         and not np.isnan(r["final_mean_fitness"])]
            mff_won = float(np.mean(won_fits)) if won_fits else float("nan")
            mff_lost = float(np.mean(lost_fits)) if lost_fits else float("nan")
            ffs = [r["final_frequency"] for r in mode1_rows]
            row = {
                "implementation": impl,
                "n1_init": n1,
                "n6_init": K - n1,
                "n_replicates": n_reps,
                "p_mode1_wins_gt_95": p_wins,
                "p_mode1_founder_extinct_within_50": p_founder_ext,
                "median_crossover_gen": median_cross,
                "mean_final_fitness_mode1_winning": mff_won,
                "mean_final_fitness_mode1_losing": mff_lost,
                "mean_final_freq_mode1": float(np.mean(ffs)),
                "sd_final_freq_mode1": float(np.std(ffs)),
                "ci_lo_p_wins": ci_lo,
                "ci_hi_p_wins": ci_hi,
            }
            if include_h:
                row["h"] = h
            w.writerow(row)


def wilsonCi(k, n, z=1.96):
    """Wilson score 95% CI for binomial proportion."""
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return float(max(0.0, centre - half)), float(min(1.0, centre + half))


# ============================================================================
# plotting (called after sweep finishes)
# ============================================================================
def plotFounderBoundary(rows, path):
    """P(Mode 1 wins) vs N_1 for Implementation A with 95% CI shading."""
    a_rows = [r for r in rows if r["implementation"] == "A"]
    cells = {}
    for r in a_rows:
        cells.setdefault(r["n1_init"], []).append(r)
    n1_sorted = sorted(cells.keys())
    p_wins, ci_lo, ci_hi = [], [], []
    for n1 in n1_sorted:
        mode1_rs = [r for r in cells[n1] if r["mode_label"] == "Mode1"]
        n = len(mode1_rs)
        k = sum(1 for r in mode1_rs if r["final_frequency"] > 0.95)
        p_wins.append(k / n if n > 0 else 0.0)
        lo, hi = wilsonCi(k, n)
        ci_lo.append(lo)
        ci_hi.append(hi)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(n1_sorted, p_wins, "o-", color="#2563eb", linewidth=2.0,
            label="Implementation A")
    ax.fill_between(n1_sorted, ci_lo, ci_hi, color="#2563eb", alpha=0.2)
    ax.axhline(0.5, color="grey", linestyle=":", alpha=0.6, label="P=0.5")
    ax.axhline(0.95, color="grey", linestyle="--", alpha=0.6, label="P=0.95")
    ax.set_xlabel("Mode 1 initial count N_1 (out of K=400)")
    ax.set_ylabel("P(Mode 1 reaches >95% by gen 1000)")
    ax.set_title("Test H2 -- Founder-loss boundary (Implementation A, 30 reps/cell, 95% CI)")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotImplementationComparison(rows, path):
    """A vs B side-by-side."""
    fig, ax = plt.subplots(figsize=(8, 5))
    color_map = {"A": "#2563eb", "B": "#ea580c"}
    for impl in ("A", "B"):
        sub = [r for r in rows if r["implementation"] == impl]
        cells = {}
        for r in sub:
            cells.setdefault(r["n1_init"], []).append(r)
        n1_sorted = sorted(cells.keys())
        p_wins, ci_lo, ci_hi = [], [], []
        for n1 in n1_sorted:
            mode1_rs = [r for r in cells[n1] if r["mode_label"] == "Mode1"]
            n = len(mode1_rs)
            k = sum(1 for r in mode1_rs if r["final_frequency"] > 0.95)
            p_wins.append(k / n if n > 0 else 0.0)
            lo, hi = wilsonCi(k, n)
            ci_lo.append(lo)
            ci_hi.append(hi)
        ax.plot(n1_sorted, p_wins, "o-", color=color_map[impl], linewidth=2.0,
                label=f"Implementation {impl}")
        ax.fill_between(n1_sorted, ci_lo, ci_hi, color=color_map[impl], alpha=0.18)
    ax.axhline(0.5, color="grey", linestyle=":", alpha=0.6)
    ax.axhline(0.95, color="grey", linestyle="--", alpha=0.6)
    ax.set_xlabel("Mode 1 initial count N_1 (out of K=400)")
    ax.set_ylabel("P(Mode 1 reaches >95% by gen 1000)")
    ax.set_title("Test H2 -- A vs B founder-loss boundary (30 reps/cell, 95% CI)")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotGradedCopyability(rows, path):
    """P(Mode 1 wins) vs N_1, one curve per h (Implementation C)."""
    c_rows = [r for r in rows if r["implementation"] == "C"]
    fig, ax = plt.subplots(figsize=(9, 6))
    cmap = plt.get_cmap("viridis")
    h_values = sorted({r["h"] for r in c_rows})
    for i, h in enumerate(h_values):
        sub = [r for r in c_rows if r["h"] == h]
        cells = {}
        for r in sub:
            cells.setdefault(r["n1_init"], []).append(r)
        n1_sorted = sorted(cells.keys())
        p_wins, ci_lo, ci_hi = [], [], []
        for n1 in n1_sorted:
            mode1_rs = [r for r in cells[n1] if r["mode_label"] == "Mode1"]
            n = len(mode1_rs)
            k = sum(1 for r in mode1_rs if r["final_frequency"] > 0.95)
            p_wins.append(k / n if n > 0 else 0.0)
            lo, hi = wilsonCi(k, n)
            ci_lo.append(lo)
            ci_hi.append(hi)
        color = cmap(i / max(1, len(h_values) - 1))
        ax.plot(n1_sorted, p_wins, "o-", color=color, linewidth=2.0,
                label=f"h={h:.2f}")
        ax.fill_between(n1_sorted, ci_lo, ci_hi, color=color, alpha=0.12)
    ax.axhline(0.5, color="grey", linestyle=":", alpha=0.6)
    ax.axhline(0.95, color="grey", linestyle="--", alpha=0.6)
    ax.set_xlabel("Mode 1 initial count N_1 (out of K=400)")
    ax.set_ylabel("P(Mode 1 reaches >95% by gen 1000)")
    ax.set_title("Test H2 -- Graded copyability (Implementation C, 30 reps/cell)")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ============================================================================
# predictions evaluator (writes test_h2_predictions_v1.md after sweep)
# ============================================================================
def evaluateAndWritePredictions(rows, path):
    """Evaluate P_H2_1..P_H2_5 and write the markdown report."""
    a_rows = [r for r in rows if r["implementation"] == "A"]
    b_rows = [r for r in rows if r["implementation"] == "B"]
    c_rows = [r for r in rows if r["implementation"] == "C"]

    # ---- helper: per-(impl, h, N_1) cell P(Mode 1 wins >95%)
    def pWinsByCell(sub_rows, impl, h=None):
        cells = {}
        for r in sub_rows:
            if r["implementation"] != impl:
                continue
            if h is not None and r["h"] != h:
                continue
            cells.setdefault(r["n1_init"], []).append(r)
        out = {}
        for n1, rs in cells.items():
            mode1_rs = [r for r in rs if r["mode_label"] == "Mode1"]
            n = len(mode1_rs)
            k = sum(1 for r in mode1_rs if r["final_frequency"] > 0.95)
            out[n1] = (k / n if n > 0 else 0.0, n, k)
        return out

    # ---- P_H2_1: monotonic + thresholds
    a_p = pWinsByCell(a_rows, "A")
    n1_sorted = sorted(a_p.keys())
    p_seq = [a_p[n1][0] for n1 in n1_sorted]
    monotonic = all(p_seq[i] <= p_seq[i + 1] + 0.10 for i in range(len(p_seq) - 1))
    p_at_20 = a_p.get(20, (None,))[0]
    p_at_80 = a_p.get(80, (None,))[0]
    p_h2_1_pass = (monotonic and (p_at_20 is None or p_at_20 < 0.5)
                   and (p_at_80 is None or p_at_80 > 0.95))
    p_h2_1_emp = (f"P(Mode 1 wins) sequence by N_1: " +
                  ", ".join(f"{n1}:{a_p[n1][0]:.2f}" for n1 in n1_sorted) +
                  f"; monotonic_within_noise={monotonic}; "
                  f"p@N_1=20={p_at_20}; p@N_1=80={p_at_80}")
    p_h2_1_verdict = "CONFIRMED" if p_h2_1_pass else "FALSIFIED"

    # ---- P_H2_2: A robustness (this is just a re-statement that A's boundary
    # exists and is well-defined). We report the threshold N_1 at which P_wins
    # crosses 0.5 and 0.95.
    cross_50_a = findCrossover(a_p, 0.5)
    cross_95_a = findCrossover(a_p, 0.95)
    p_h2_2_emp = (f"Implementation A boundary: N_1 at P_wins=0.5 is ~{cross_50_a}, "
                  f"at P_wins=0.95 is ~{cross_95_a}")
    # P_H2_2 is "CONFIRMED" if a well-defined boundary exists (interpolated cross)
    p_h2_2_verdict = "CONFIRMED" if cross_50_a is not None and cross_95_a is not None else "INCONCLUSIVE"

    # ---- P_H2_3: B matches A
    b_p = pWinsByCell(b_rows, "B")
    common_n1 = sorted(set(a_p.keys()) & set(b_p.keys()))
    diffs = [abs(a_p[n1][0] - b_p[n1][0]) for n1 in common_n1]
    max_diff = max(diffs) if diffs else float("nan")
    # threshold for "match": all per-cell diffs < 0.2 (within Wilson 95% CI for 30 reps)
    p_h2_3_pass = max_diff < 0.2 if diffs else False
    p_h2_3_emp = (f"max per-cell |P_A - P_B| = {max_diff:.3f} across {len(common_n1)} shared N_1 cells; "
                  f"per-cell A vs B: " +
                  ", ".join(f"N_1={n1}:A={a_p[n1][0]:.2f}/B={b_p[n1][0]:.2f}" for n1 in common_n1))
    p_h2_3_verdict = "CONFIRMED" if p_h2_3_pass else "FALSIFIED"

    # ---- P_H2_4: C with h=0 ~ A, h=1 ~ Mode 1, h=0.5 still loses
    c_h_values = sorted({r["h"] for r in c_rows})
    c_results = {h: pWinsByCell(c_rows, "C", h=h) for h in c_h_values}
    # at h=0.5, check Mode 1 still wins at low N_1 like A does (or worse)
    h05 = 0.5
    if h05 in c_results:
        p05 = c_results[h05]
        cross_50_h05 = findCrossover(p05, 0.5)
        # at h=1.0, behaves like Mode 1 vs Mode 1 -> P_wins at N_1=200 should be ~0.5
        p_at_200_h1 = c_results.get(1.0, {}).get(200, (None,))[0]
        p_h2_4_emp = (f"crossover N_1 at h=0.5 is ~{cross_50_h05}; "
                      f"P(Mode 1 wins) at N_1=200, h=1.0 is {p_at_200_h1}")
    else:
        cross_50_h05 = None
        p_h2_4_emp = "h=0.5 cell missing"
    # framework's prediction: at h=0.5, crossover is closer to A than to N_1=200
    p_h2_4_pass = (cross_50_h05 is not None and cross_50_a is not None
                   and abs(cross_50_h05 - cross_50_a) < (200 - cross_50_a) * 0.6)
    p_h2_4_verdict = "CONFIRMED" if p_h2_4_pass else "REFINED"

    # ---- P_H2_5: graded crossover shifts smoothly with h
    crossovers_by_h = {}
    for h in c_h_values:
        crossovers_by_h[h] = findCrossover(c_results[h], 0.5)
    cross_seq = [(h, crossovers_by_h[h]) for h in sorted(crossovers_by_h.keys())]
    p_h2_5_emp = ("crossover N_1 (where P_wins=0.5) by h: " +
                  ", ".join(f"h={h:.2f}:N_1~{c}" for h, c in cross_seq))
    # check shifting right with h
    cross_vals = [c for _, c in cross_seq if c is not None]
    if len(cross_vals) >= 2:
        shifts_right = all(cross_vals[i] <= cross_vals[i + 1] + 30
                           for i in range(len(cross_vals) - 1))
        p_h2_5_pass = shifts_right
    else:
        p_h2_5_pass = False
    p_h2_5_verdict = "CONFIRMED" if p_h2_5_pass else "REFINED"

    # ---- write markdown
    lines = []
    lines.append("# Test H2 -- Pre-registered prediction results (v1)")
    lines.append("")
    lines.append("Pre-registered in `code/test_h2_competition_sweep.py` header BEFORE simulation runs.")
    lines.append("")
    lines.append("## Predictions table")
    lines.append("")
    lines.append("| ID      | Verdict      | Empirical observation |")
    lines.append("|---------|--------------|------------------------|")
    for pid, verdict, emp in [
        ("P_H2_1", p_h2_1_verdict, p_h2_1_emp),
        ("P_H2_2", p_h2_2_verdict, p_h2_2_emp),
        ("P_H2_3", p_h2_3_verdict, p_h2_3_emp),
        ("P_H2_4", p_h2_4_verdict, p_h2_4_emp),
        ("P_H2_5", p_h2_5_verdict, p_h2_5_emp),
    ]:
        emp_safe = emp.replace("|", "\\|")
        lines.append(f"| {pid} | **{verdict}** | {emp_safe} |")
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append(f"- Sweep totals: Implementation A = {len(STEP1_N1_VALUES)} N_1 cells x {N_REPLICATES} reps; "
                 f"Implementation B = {len(STEP2_N1_VALUES)} cells x {N_REPLICATES} reps; "
                 f"Implementation C = {len(STEP3_H_VALUES)} h x {len(STEP3_N1_VALUES)} N_1 cells x {N_REPLICATES} reps.")
    lines.append(f"- K = {K}, L_TARGET = {L_TARGET}, N_GEN = {N_GEN}, BETA = {BETA}, MU_MODE1 = {MU_MODE1}, EPS_NOISE_6 = {EPS_NOISE_6}")
    lines.append("- 95% CI from Wilson score on the per-cell binomial.")
    lines.append("- Seeding: seed = 42 + cell_idx*100 + rep_idx (unique per (cell, rep) globally).")
    Path(path).write_text("\n".join(lines))


def findCrossover(p_dict, threshold):
    """Find N_1 where P_wins crosses `threshold` (linear interp between cells).

    p_dict: {n1: (p_wins, n, k)}
    Returns interpolated N_1 or None if no crossing.
    """
    n1_sorted = sorted(p_dict.keys())
    p_seq = [p_dict[n1][0] for n1 in n1_sorted]
    for i in range(len(n1_sorted) - 1):
        p0, p1 = p_seq[i], p_seq[i + 1]
        if (p0 <= threshold <= p1) or (p1 <= threshold <= p0):
            if p1 == p0:
                return n1_sorted[i]
            frac = (threshold - p0) / (p1 - p0)
            return float(n1_sorted[i] + frac * (n1_sorted[i + 1] - n1_sorted[i]))
    return None


# ============================================================================
# signal handler for clean shutdown
# ============================================================================
def installSignalHandlers(state):
    def handler(signum, frame):
        state["shutdown"] = True
        writeFinalProgress(state, status=f"signal_{signum}_received")
        sys.exit(0)
    signal.signal(signal.SIGTERM, handler)
    signal.signal(signal.SIGINT, handler)


# ============================================================================
# main
# ============================================================================
def main():
    np.random.seed(42)

    parser = argparse.ArgumentParser()
    parser.add_argument("--progress-file", required=True,
                        help="path to write progress lines")
    parser.add_argument("--completed-csv", required=True,
                        help="path to append every finished sim row")
    parser.add_argument("--smoke-test", action="store_true",
                        help="tiny test run for verification (Implementation A only, "
                             "2 N_1 values x 2 reps x 50 generations)")
    args = parser.parse_args()

    progress_path = Path(args.progress_file)
    completed_path = Path(args.completed_csv)
    progress_path.parent.mkdir(parents=True, exist_ok=True)
    completed_path.parent.mkdir(parents=True, exist_ok=True)

    # truncate progress.txt at start
    progress_path.write_text(f"# === test_h2 run started at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")

    if args.smoke_test:
        # tiny config: Implementation A only, 2 cells x 2 reps x 50 gens
        n_replicates = 2
        n_gen = 50
        smoke_cells = [("A", 0.0, 40), ("A", 0.0, 200)]
        cells = smoke_cells
        # truncate completed.csv for smoke
        if completed_path.exists():
            completed_path.unlink()
    else:
        n_replicates = N_REPLICATES
        n_gen = N_GEN
        cells = buildCellList()
        # if completed.csv exists, archive it first to avoid contamination
        if completed_path.exists():
            backup = completed_path.with_suffix(".csv.bak")
            completed_path.rename(backup)
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] archived prior completed.csv to {backup}\n")

    n_total = len(cells) * n_replicates
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
                f"{len(cells)} cells x {n_replicates} reps = {n_total} sims; "
                f"n_gen={n_gen}; smoke_test={args.smoke_test}\n")

    # run all cells
    for cell_idx, (impl, h, n1) in enumerate(cells):
        if state["shutdown"]:
            break
        state["current_step"] = f"impl={impl},h={h:.2f},N_1={n1}"
        runCell(impl, h, n1, cell_idx, n_replicates, n_gen, completed_path, state)
        # write a progress line at end of every cell as well (in addition to every-50 rule)
        writeProgressLine(state)

    # final progress line
    writeFinalProgress(state, status="completed")

    # post-process: build the per-step CSVs, plots, and predictions markdown
    # (only for the full run; skip for smoke)
    if not args.smoke_test:
        try:
            all_rows = loadCompletedRows(completed_path)
            writeStepCsvs(all_rows)
            plotFounderBoundary(all_rows, FIGURES_DIR / "test_h2_founder_boundary.png")
            plotImplementationComparison(all_rows, FIGURES_DIR / "test_h2_implementation_comparison.png")
            plotGradedCopyability(all_rows, FIGURES_DIR / "test_h2_graded_copyability.png")
            evaluateAndWritePredictions(all_rows, RESULTS_DIR / "test_h2_predictions_v1.md")
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] post-process complete: "
                        f"3 step CSVs, 3 figures, predictions_v1.md written\n")
        except Exception as e:
            with open(progress_path, "a") as f:
                f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] post-process FAILED: {e}\n")
            raise


if __name__ == "__main__":
    main()
