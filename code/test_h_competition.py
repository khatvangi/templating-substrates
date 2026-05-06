"""
Test H -- Origin-of-life Mode 1 vs Mode 6 competition (and Mode 5 head-to-heads).

Tests the framework's strongest causal claim, dispatched from the published draft:
section 4 (line 127) of templating_substrates_draft_v2.md asserts that "an
origin-of-life simulation with multiple modes equally available should not
develop Mode 6-based inheritance." Test H is the literal pre-registered
execution of that prediction. The framework's claim is that Mode 6's lack of a
complementarity-based copy mechanism causes its competitive extinction within a
bounded number of generations *regardless of initial abundance*. Test D
established each mode's solo fitness ceiling; Test H couples the populations in
a shared substrate pool under shared selection and watches who wins.

================================================================================
Pre-registered predictions (DEFINED BEFORE RUN, see dispatch3.md)
================================================================================

P_H1: Starting with equal-population Mode 1 and Mode 6 (200 agents each),
      Mode 6 frequency declines monotonically and reaches < 5% of total
      population within 200 generations under selection toward a fixed
      random target sequence of length 32.

P_H2: The crossover generation (Mode 1 frequency exceeds 95%) is bounded
      above by 500 generations, and the crossover dynamics are dominated
      by Mode 1's adaptation rate, not by Mode 6's slower extinction.

P_H3: When Mode 6 starts with 90% of initial population (Mode 1 at 10%),
      Mode 1 still climbs to > 95% by generation 1000. The framework's
      claim is that initial abundance does not matter -- copyability does.
      P_H3 is the strongest version of the framework's claim. Failure here
      would show that initial conditions matter more than the framework
      admits.

P_H4: When Mode 1 is replaced with Mode 5 (length-bounded inheritance),
      the competition is closer: Mode 5 dominates Mode 6 but plateaus at
      the Mode 5 ceiling (~ 0.41 fitness), and may not reach > 95%
      depending on Mode 6's initial abundance. Mode 5 vs Mode 1 still
      favors Mode 1 because Mode 1's ceiling is higher.

================================================================================
Scenarios (each: 10 replicates, 1000 generations, K = 400)
================================================================================

A: Mode 1 (200) vs Mode 6 (200)             -- equal start, P_H1 / P_H2
B: Mode 1 ( 40) vs Mode 6 (360)             -- Mode 1 starts at 10%, P_H3
C: Mode 5 (200) vs Mode 6 (200)             -- Mode 5 head-to-head, P_H4 (low-info)
D: Mode 1 (200) vs Mode 5 (200)             -- inheritance-capable head-to-head, P_H4 (high-info)

================================================================================
Mode 6 implementation choice (matters -- the framework's claim hinges on it)
================================================================================

This script uses Implementation A (binary copyability) per dispatch3.md
Step 1. Each Mode 6 agent draws phenotype = template + i.i.d. noise where the
TEMPLATE is fixed at the start of the run for the whole Mode 6 population (a
single 2D-surface-pattern that all Mode 6 agents produce noisy realizations
of, with NO inheritance of adaptation). Mutations do not propagate because
there is no copy mechanism. This matches the dispatch's gloss: "Mode 6
produces consistent patterns the population cannot evolve away from."

Implementations B (per-lineage fixed template) and C (heritability < 1) from
the dispatch §Step 4 are not implemented in this v1; deferred to a sensitivity
follow-up if A's results warrant it.

================================================================================
Reproducibility
================================================================================

Module-level seed: np.random.seed(42) and np.random.default_rng(42).
Per-replicate seeds: replicate r in {0..9} uses seed 42 + r, i.e. {42..51}.

================================================================================
Parameters (matched to test_d_v2)
================================================================================

ALPHABET    = 4
L_TARGET    = 32
K           = 400  (total population cap, shared substrate)
N_GEN       = 1000
BETA        = 10.0 (selection strength)
MU_MODE1    = 0.01  per-position
MU_MODE5    = 0.02  per-position (matches test_d_v2)
N_MODULES_5 = 8
EPS_NOISE_5 = 0.05  (Mode 5 phenotype noise on templated positions)
EPS_NOISE_6 = 0.05  (Mode 6 phenotype noise around the fixed template)

================================================================================
Anti-DRY discipline
================================================================================

Per the project's CLAUDE.md (rule 4), reproduction / mutation / selection
logic is duplicated into this script rather than imported from
test_d_v2_population_dynamics.py. If you fix a bug in one, fix it
intentionally in the other.
"""

import csv
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
# parameters (must match test_d_v2)
# ============================================================================
ALPHABET = 4
L_TARGET = 32
K = 400
N_GEN = 1000
BETA = 10.0

MU_MODE1 = 0.01
MU_MODE5 = 0.02
N_MODULES_5 = 8
EPS_NOISE_5 = 0.05
EPS_NOISE_6 = 0.05

N_REPLICATES = 10

# fixed target across the whole run (same as test_d_v2 idiom: seed 2026)
TARGET = np.random.default_rng(2026).integers(0, ALPHABET, size=L_TARGET, dtype=np.int8)


# ============================================================================
# fitness (vectorized over agent pool)
# ============================================================================
def findFitness(phenotypes):
    """phenotypes shape (N, L_TARGET) -> fitness of shape (N,) in [0,1]."""
    matches = (phenotypes == TARGET[None, :]).sum(axis=1)
    return matches.astype(np.float64) / L_TARGET


# ============================================================================
# Mode 1: 1D template, faithful copy with point mutation rate MU_MODE1
# duplicated from test_d_v2_population_dynamics.py per anti-DRY discipline
# ============================================================================
class Mode1Subpop:
    label = "Mode1"

    def __init__(self, n_agents, rng):
        # genotype = phenotype (length L_TARGET)
        self.rng = rng
        self.templates = rng.integers(0, ALPHABET, size=(n_agents, L_TARGET), dtype=np.int8)

    @property
    def n(self):
        return self.templates.shape[0]

    def findPhenotypes(self):
        # mode 1 phenotype is the template itself (no noise per test_d_v2)
        return self.templates

    def reproduceFromIndices(self, parent_indices):
        # copy parent template with per-position mutation rate MU_MODE1
        new_templates = self.templates[parent_indices].copy()
        flips = self.rng.random(new_templates.shape) < MU_MODE1
        replacements = self.rng.integers(0, ALPHABET, size=new_templates.shape, dtype=np.int8)
        new_templates = np.where(flips, replacements, new_templates)
        self.templates = new_templates


# ============================================================================
# Mode 5: modular conveyor with N_MODULES_5 modules
# duplicated from test_d_v2_population_dynamics.py per anti-DRY discipline
# ============================================================================
class Mode5Subpop:
    label = f"Mode5_N{N_MODULES_5}"

    def __init__(self, n_agents, rng):
        self.rng = rng
        # genotype = N_MODULES_5 module preferences
        self.module_prefs = rng.integers(0, ALPHABET, size=(n_agents, N_MODULES_5), dtype=np.int8)

    @property
    def n(self):
        return self.module_prefs.shape[0]

    def findPhenotypes(self):
        # first N_MODULES_5 positions are templated (with EPS_NOISE_5);
        # remaining tail positions are uniform random per draw (untemplated)
        n_agents = self.module_prefs.shape[0]
        out = np.empty((n_agents, L_TARGET), dtype=np.int8)
        L_temp = min(N_MODULES_5, L_TARGET)
        intended = self.module_prefs[:, :L_temp]
        noise = self.rng.random(intended.shape) < EPS_NOISE_5
        random_subs = self.rng.integers(0, ALPHABET, size=intended.shape, dtype=np.int8)
        out[:, :L_temp] = np.where(noise, random_subs, intended)
        if L_TARGET > N_MODULES_5:
            out[:, N_MODULES_5:] = self.rng.integers(
                0, ALPHABET, size=(n_agents, L_TARGET - N_MODULES_5), dtype=np.int8
            )
        return out

    def reproduceFromIndices(self, parent_indices):
        # copy parent module prefs with mutation rate MU_MODE5
        new_prefs = self.module_prefs[parent_indices].copy()
        flips = self.rng.random(new_prefs.shape) < MU_MODE5
        replacements = self.rng.integers(0, ALPHABET, size=new_prefs.shape, dtype=np.int8)
        new_prefs = np.where(flips, replacements, new_prefs)
        self.module_prefs = new_prefs


# ============================================================================
# Mode 6 -- Implementation A (binary copyability):
#   * Mode 6 has a single FIXED template (the 2D-surface pattern) shared by
#     the entire Mode 6 population for the entire run.
#   * Phenotype of each agent each generation = template + i.i.d. noise per
#     position (noise rate EPS_NOISE_6).
#   * Reproduction does not propagate any inherited variation; we keep a
#     trivial per-agent "genotype" for bookkeeping, but it is irrelevant to
#     the phenotype distribution. Operationally: when an agent is selected
#     to reproduce, the offspring just gets a fresh draw from the same
#     fixed-template noise process. There is no copy mechanism, so
#     adaptation cannot accumulate.
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
        # i.i.d. draws around the fixed template
        n_agents = self._n
        if n_agents == 0:
            return np.empty((0, L_TARGET), dtype=np.int8)
        intended = np.broadcast_to(self.fixed_template[None, :], (n_agents, L_TARGET))
        noise = self.rng.random((n_agents, L_TARGET)) < EPS_NOISE_6
        random_subs = self.rng.integers(0, ALPHABET, size=(n_agents, L_TARGET), dtype=np.int8)
        return np.where(noise, random_subs, intended).astype(np.int8)

    def reproduceFromIndices(self, parent_indices):
        # no genotype, no copying: the only state is the count.
        # whichever parent indices are picked, the next generation's agents
        # will draw fresh phenotypes from the same fixed template.
        self._n = len(parent_indices)


# ============================================================================
# competition driver: two subpopulations sharing K slots
# ============================================================================
def runCompetition(make_subpop_a, make_subpop_b, n_a_init, n_b_init, n_gen, rng):
    """
    Runs a coupled competition for n_gen generations.

    make_subpop_a / make_subpop_b: callables (n_agents, rng) -> Subpop
    n_a_init / n_b_init: initial counts (must sum to K)
    rng: numpy Generator seeded for the replicate

    Returns: list of dicts with per-generation per-mode stats.
    """
    assert n_a_init + n_b_init == K, f"initial counts must sum to K={K}"

    # we use distinct rng substreams per subpop to keep within-subpop
    # mutations independent of the cross-population selection draw.
    rng_a = np.random.default_rng(rng.integers(0, 2**31 - 1))
    rng_b = np.random.default_rng(rng.integers(0, 2**31 - 1))

    sub_a = make_subpop_a(n_a_init, rng_a)
    sub_b = make_subpop_b(n_b_init, rng_b)

    history = []

    for g in range(n_gen):
        # 1. each subpop produces its phenotypes for this generation
        phen_a = sub_a.findPhenotypes() if sub_a.n > 0 else np.empty((0, L_TARGET), dtype=np.int8)
        phen_b = sub_b.findPhenotypes() if sub_b.n > 0 else np.empty((0, L_TARGET), dtype=np.int8)

        fit_a = findFitness(phen_a) if sub_a.n > 0 else np.empty(0, dtype=np.float64)
        fit_b = findFitness(phen_b) if sub_b.n > 0 else np.empty(0, dtype=np.float64)

        # 2. record stats BEFORE selection (this is the state at gen g)
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

        # 3. shared selection: stack fitnesses, weight by exp(BETA * fit), draw K parents
        all_fit = np.concatenate([fit_a, fit_b])
        if all_fit.size == 0:
            break
        # numerical-stable softmax-like weighting
        w = np.exp(BETA * (all_fit - all_fit.max()))
        w = w / w.sum()
        parent_indices = rng.choice(n_total, size=K, replace=True, p=w)

        # 4. partition selected parents back into the two subpops
        is_a = parent_indices < n_a
        a_parents = parent_indices[is_a]               # indices into sub_a
        b_parents = parent_indices[~is_a] - n_a        # indices into sub_b (offset back)

        # 5. reproduce within each subpop using its own selected parents
        if a_parents.size > 0:
            sub_a.reproduceFromIndices(a_parents)
        else:
            # extinct: zero out
            sub_a = make_subpop_a(0, rng_a)
        if b_parents.size > 0:
            sub_b.reproduceFromIndices(b_parents)
        else:
            sub_b = make_subpop_b(0, rng_b)

    return history


# ============================================================================
# scenario factories
# ============================================================================
def makeMode1(n_agents, rng):
    return Mode1Subpop(n_agents, rng)


def makeMode5(n_agents, rng):
    return Mode5Subpop(n_agents, rng)


def makeMode6_factory(fixed_template):
    """returns a closure (n, rng) -> Mode6SubpopA bound to one fixed template."""
    def maker(n_agents, rng):
        return Mode6SubpopA(n_agents, rng, fixed_template)
    return maker


# ============================================================================
# replicate runner per scenario
# ============================================================================
def runScenarioReplicates(scenario_name, make_a, make_b_factory_or_callable, n_a_init, n_b_init,
                          uses_mode6, progress_writer, n_replicates=N_REPLICATES, n_gen=N_GEN):
    """
    Runs n_replicates of one scenario.

    make_b_factory_or_callable: if uses_mode6, this is a *factory* that needs
    a fixed template per replicate; otherwise it is a regular maker.
    """
    all_rows = []
    for r in range(n_replicates):
        seed = 42 + r
        rng = np.random.default_rng(seed)
        # build mode 6's fixed template per replicate (so reps differ) using a
        # derived stream from the master rng so we keep determinism with seed.
        if uses_mode6:
            tmpl_rng = np.random.default_rng(seed * 7919 + 1)  # arbitrary mixing
            fixed_template = tmpl_rng.integers(0, ALPHABET, size=L_TARGET, dtype=np.int8)
            make_b = make_b_factory_or_callable(fixed_template)
        else:
            make_b = make_b_factory_or_callable

        t0 = time.time()
        hist = runCompetition(make_a, make_b, n_a_init, n_b_init, n_gen, rng)
        dt = time.time() - t0

        # tag rows with replicate + scenario
        for row in hist:
            row["replicate"] = r
            row["seed"] = seed
            row["scenario"] = scenario_name
            all_rows.append(row)

        line = (f"[{time.strftime('%H:%M:%S')}] scenario={scenario_name} rep={r} "
                f"seed={seed} dt={dt:.1f}s "
                f"final_freqs=({hist[-2]['frequency']:.3f}, {hist[-1]['frequency']:.3f})")
        progress_writer.write(line + "\n")
        progress_writer.flush()
        print(line)

    return all_rows


# ============================================================================
# helpers: per-scenario summary statistics
# ============================================================================
def summarizeScenario(rows, scenario_name, mode_a_label, mode_b_label):
    """
    Returns one dict per (replicate, mode) summarising the run.

    Fields:
      scenario, replicate, mode_label, final_count, final_frequency,
      final_mean_fitness, max_mean_fitness, generation_of_extinction (<5%),
      generation_of_dominance (>95%)
    """
    out = []
    by_rep_mode = {}
    for row in rows:
        if row["scenario"] != scenario_name:
            continue
        key = (row["replicate"], row["mode_label"])
        by_rep_mode.setdefault(key, []).append(row)

    for (rep, mode), rows_rm in by_rep_mode.items():
        rows_rm.sort(key=lambda r: r["generation"])
        gens = [r["generation"] for r in rows_rm]
        freqs = [r["frequency"] for r in rows_rm]
        fits = [r["mean_fitness"] for r in rows_rm]

        # generation of extinction = first gen at which freq < 0.05
        ext_gen = None
        for g, f in zip(gens, freqs):
            if f < 0.05:
                ext_gen = g
                break
        # generation of dominance = first gen at which freq > 0.95
        dom_gen = None
        for g, f in zip(gens, freqs):
            if f > 0.95:
                dom_gen = g
                break

        # max mean fitness ignoring NaNs (subpop may have gone extinct)
        finite_fits = [x for x in fits if not (x is None or np.isnan(x))]
        out.append({
            "scenario": scenario_name,
            "replicate": rep,
            "mode_label": mode,
            "final_count": rows_rm[-1]["count"],
            "final_frequency": rows_rm[-1]["frequency"],
            "final_mean_fitness": rows_rm[-1]["mean_fitness"],
            "max_mean_fitness": max(finite_fits) if finite_fits else float("nan"),
            "generation_of_extinction_lt_5pct": ext_gen if ext_gen is not None else -1,
            "generation_of_dominance_gt_95pct": dom_gen if dom_gen is not None else -1,
        })
    return out


# ============================================================================
# CSV writers
# ============================================================================
def writeScenarioCsv(rows, scenario_name, path):
    fields = ["scenario", "replicate", "seed", "generation", "mode_label",
              "count", "frequency", "mean_fitness", "max_fitness",
              "total_population"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            if row["scenario"] != scenario_name:
                continue
            w.writerow({k: row.get(k, "") for k in fields})


def writeSummaryCsv(summary_rows, path):
    fields = ["scenario", "replicate", "mode_label",
              "final_count", "final_frequency", "final_mean_fitness",
              "max_mean_fitness",
              "generation_of_extinction_lt_5pct",
              "generation_of_dominance_gt_95pct"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in summary_rows:
            w.writerow({k: row[k] for k in fields})


# ============================================================================
# plotting
# ============================================================================
def plotScenarioFreqVsGen(rows, scenario_name, path):
    """one panel: mean +/- SD of frequency vs generation per mode."""
    by_mode = {}
    for row in rows:
        if row["scenario"] != scenario_name:
            continue
        by_mode.setdefault(row["mode_label"], {}).setdefault(row["generation"], []).append(row["frequency"])

    fig, ax = plt.subplots(figsize=(8, 5))
    color_map = {"Mode1": "#2563eb", f"Mode5_N{N_MODULES_5}": "#ea580c",
                 "Mode6_implA": "#0891b2"}
    for mode, gen_dict in by_mode.items():
        gens = sorted(gen_dict.keys())
        means = np.array([np.mean(gen_dict[g]) for g in gens])
        sds = np.array([np.std(gen_dict[g]) for g in gens])
        color = color_map.get(mode, "black")
        ax.plot(gens, means, "-", color=color, linewidth=2.0, label=mode)
        ax.fill_between(gens, means - sds, means + sds, color=color, alpha=0.2)
    ax.axhline(0.95, color="grey", linestyle=":", alpha=0.5, label="95%")
    ax.axhline(0.05, color="grey", linestyle=":", alpha=0.5)
    ax.set_xlabel("generation")
    ax.set_ylabel("population frequency")
    ax.set_title(f"Test H scenario {scenario_name} -- frequency vs generation (mean +/- SD across reps)")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="center right")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plotAllScenarios(rows, path):
    """small multiples: 2x2 grid, one panel per scenario."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharey=True)
    color_map = {"Mode1": "#2563eb", f"Mode5_N{N_MODULES_5}": "#ea580c",
                 "Mode6_implA": "#0891b2"}
    scenario_layout = [("A", axes[0][0]), ("B", axes[0][1]),
                       ("C", axes[1][0]), ("D", axes[1][1])]
    for scen, ax in scenario_layout:
        by_mode = {}
        for row in rows:
            if row["scenario"] != scen:
                continue
            by_mode.setdefault(row["mode_label"], {}).setdefault(row["generation"], []).append(row["frequency"])
        for mode, gen_dict in by_mode.items():
            gens = sorted(gen_dict.keys())
            means = np.array([np.mean(gen_dict[g]) for g in gens])
            sds = np.array([np.std(gen_dict[g]) for g in gens])
            color = color_map.get(mode, "black")
            ax.plot(gens, means, "-", color=color, linewidth=1.8, label=mode)
            ax.fill_between(gens, means - sds, means + sds, color=color, alpha=0.2)
        ax.axhline(0.95, color="grey", linestyle=":", alpha=0.4)
        ax.axhline(0.05, color="grey", linestyle=":", alpha=0.4)
        ax.set_title(f"scenario {scen}")
        ax.set_xlabel("generation")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="center right")
    axes[0][0].set_ylabel("population frequency")
    axes[1][0].set_ylabel("population frequency")
    fig.suptitle("Test H -- four competition scenarios (mean +/- SD across 10 replicates)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ============================================================================
# prediction evaluator
# ============================================================================
def evaluatePredictions(summary_rows):
    """
    Returns a dict keyed by prediction id with:
      - statement (string)
      - empirical_value (string)
      - verdict (CONFIRMED / FALSIFIED)
      - notes (string)
    """
    def repsForScenario(scen, mode):
        return [r for r in summary_rows if r["scenario"] == scen and r["mode_label"] == mode]

    res = {}

    # ---- P_H1: scenario A, Mode 6 (impl A) <5% within 200 gens ----
    a_mode6 = repsForScenario("A", "Mode6_implA")
    ext_gens_h1 = [r["generation_of_extinction_lt_5pct"] for r in a_mode6]
    finite_ext_h1 = [g for g in ext_gens_h1 if g >= 0]
    n_ext_within_200 = sum(1 for g in finite_ext_h1 if g <= 200)
    n_ext_total = len(finite_ext_h1)
    if len(a_mode6) == 0:
        ph1_verdict = "INCONCLUSIVE"
        ph1_emp = "no data"
    else:
        # criterion: ALL 10 reps reach <5% by gen 200
        if n_ext_within_200 == len(a_mode6):
            ph1_verdict = "CONFIRMED"
        else:
            ph1_verdict = "FALSIFIED"
        if finite_ext_h1:
            ph1_emp = (f"{n_ext_within_200}/{len(a_mode6)} reps reach <5% by gen 200; "
                       f"median extinction gen across reps that reached <5% = "
                       f"{int(np.median(finite_ext_h1))} (range {min(finite_ext_h1)}-{max(finite_ext_h1)}); "
                       f"reps reaching <5% by end of run = {n_ext_total}/{len(a_mode6)}")
        else:
            ph1_emp = f"0/{len(a_mode6)} reps reached <5%"

    res["P_H1"] = {
        "statement": ("Scenario A: Mode 6 frequency declines to <5% within 200 generations "
                      "in all 10 replicates."),
        "empirical_value": ph1_emp,
        "verdict": ph1_verdict,
        "notes": "selection toward fixed random target len 32, K=400, 10 reps",
    }

    # ---- P_H2: scenario A, crossover (Mode 1 >95%) <= 500 gens ----
    a_mode1 = repsForScenario("A", "Mode1")
    dom_gens_h2 = [r["generation_of_dominance_gt_95pct"] for r in a_mode1]
    finite_dom_h2 = [g for g in dom_gens_h2 if g >= 0]
    if len(a_mode1) == 0:
        ph2_verdict = "INCONCLUSIVE"
        ph2_emp = "no data"
    else:
        n_dom_within_500 = sum(1 for g in finite_dom_h2 if g <= 500)
        if n_dom_within_500 == len(a_mode1):
            ph2_verdict = "CONFIRMED"
        else:
            ph2_verdict = "FALSIFIED"
        if finite_dom_h2:
            ph2_emp = (f"{n_dom_within_500}/{len(a_mode1)} reps reach Mode1>95% by gen 500; "
                       f"median crossover gen = {int(np.median(finite_dom_h2))} "
                       f"(range {min(finite_dom_h2)}-{max(finite_dom_h2)})")
        else:
            ph2_emp = f"0/{len(a_mode1)} reps reached >95%"
    res["P_H2"] = {
        "statement": "Scenario A: Mode 1 reaches >95% frequency by generation 500 in all reps.",
        "empirical_value": ph2_emp,
        "verdict": ph2_verdict,
        "notes": ("strongest version: bounded by 500 gens; weaker forms "
                  "(by 1000) reported in summary csv"),
    }

    # ---- P_H3: scenario B, Mode 1 climbs to >95% by gen 1000 ----
    b_mode1 = repsForScenario("B", "Mode1")
    dom_gens_h3 = [r["generation_of_dominance_gt_95pct"] for r in b_mode1]
    finite_dom_h3 = [g for g in dom_gens_h3 if g >= 0]
    final_freqs_h3 = [r["final_frequency"] for r in b_mode1]
    if len(b_mode1) == 0:
        ph3_verdict = "INCONCLUSIVE"
        ph3_emp = "no data"
    else:
        n_dom_h3 = sum(1 for g in finite_dom_h3 if g <= 1000)
        # criterion: ALL reps reach >95% by gen 1000
        if n_dom_h3 == len(b_mode1):
            ph3_verdict = "CONFIRMED"
        else:
            ph3_verdict = "FALSIFIED"
        ph3_emp = (f"{n_dom_h3}/{len(b_mode1)} reps reach Mode1>95% by gen 1000; "
                   f"final Mode1 frequency mean = {np.mean(final_freqs_h3):.3f}, "
                   f"median crossover gen (where reached) = "
                   f"{int(np.median(finite_dom_h3)) if finite_dom_h3 else -1}")
    res["P_H3"] = {
        "statement": ("Scenario B (Mode 1 starts at 10%): Mode 1 reaches >95% by "
                      "generation 1000 in all reps. THIS IS THE STRONGEST VERSION OF "
                      "THE FRAMEWORK'S CLAIM (initial abundance does not matter)."),
        "empirical_value": ph3_emp,
        "verdict": ph3_verdict,
        "notes": "if FALSIFIED, framework's copyability-only claim is wrong",
    }

    # ---- P_H4: scenario C (Mode 5 vs Mode 6) and scenario D (Mode 1 vs Mode 5) ----
    c_mode5 = repsForScenario("C", f"Mode5_N{N_MODULES_5}")
    c_mode6 = repsForScenario("C", "Mode6_implA")
    d_mode1 = repsForScenario("D", "Mode1")
    d_mode5 = repsForScenario("D", f"Mode5_N{N_MODULES_5}")

    if c_mode5 and c_mode6 and d_mode1 and d_mode5:
        c_m5_final = np.mean([r["final_frequency"] for r in c_mode5])
        c_m6_final = np.mean([r["final_frequency"] for r in c_mode6])
        d_m1_final = np.mean([r["final_frequency"] for r in d_mode1])
        d_m5_final = np.mean([r["final_frequency"] for r in d_mode5])
        # criterion (P_H4 has two parts):
        #   (a) C: Mode 5 dominates Mode 6 (Mode5 > Mode6 at end)
        #   (b) D: Mode 1 dominates Mode 5 (Mode1 > Mode5 at end)
        c_pass = c_m5_final > c_m6_final
        d_pass = d_m1_final > d_m5_final
        if c_pass and d_pass:
            ph4_verdict = "CONFIRMED"
        else:
            ph4_verdict = "FALSIFIED"
        ph4_emp = (f"scenario C final freqs: Mode5={c_m5_final:.3f}, Mode6={c_m6_final:.3f} "
                   f"(Mode5 dominates Mode6: {c_pass}); "
                   f"scenario D final freqs: Mode1={d_m1_final:.3f}, Mode5={d_m5_final:.3f} "
                   f"(Mode1 dominates Mode5: {d_pass})")
    else:
        ph4_verdict = "INCONCLUSIVE"
        ph4_emp = "missing scenario C or D data"

    res["P_H4"] = {
        "statement": ("Scenario C (Mode 5 vs Mode 6): Mode 5 dominates Mode 6; "
                      "Scenario D (Mode 1 vs Mode 5): Mode 1 dominates Mode 5."),
        "empirical_value": ph4_emp,
        "verdict": ph4_verdict,
        "notes": "Mode 5 may not reach 95% because of its 0.41 fitness ceiling",
    }

    return res


# ============================================================================
# crossover-generation helper for the report
# ============================================================================
def crossoverPerScenario(summary_rows):
    """
    For each scenario, return median (across reps) gen where the dominant mode
    first reaches >95%, plus mode label.
    """
    out = {}
    scenario_dominant = {
        "A": "Mode1", "B": "Mode1",
        "C": f"Mode5_N{N_MODULES_5}", "D": "Mode1",
    }
    for scen, dom_mode in scenario_dominant.items():
        rs = [r for r in summary_rows if r["scenario"] == scen and r["mode_label"] == dom_mode]
        gens = [r["generation_of_dominance_gt_95pct"] for r in rs]
        finite = [g for g in gens if g >= 0]
        n_reached = len(finite)
        n_total = len(rs)
        med = int(np.median(finite)) if finite else -1
        out[scen] = {
            "dominant_mode": dom_mode,
            "n_reps_reached_95pct": n_reached,
            "n_reps_total": n_total,
            "median_crossover_gen": med,
            "min_crossover_gen": min(finite) if finite else -1,
            "max_crossover_gen": max(finite) if finite else -1,
        }
    return out


# ============================================================================
# predictions markdown writer
# ============================================================================
def writePredictionsMd(predictions, crossovers, scenario_summary_table, path):
    lines = []
    lines.append("# Test H -- Pre-registered prediction results (v1)")
    lines.append("")
    lines.append("Pre-registered in `code/test_h_competition.py` header BEFORE simulation runs.")
    lines.append("Mode 6 implementation: A (binary copyability, fixed-template noise draw) per dispatch3.md §Step 1.")
    lines.append("")
    lines.append("## Predictions table")
    lines.append("")
    lines.append("| ID    | Verdict     | Empirical observation |")
    lines.append("|-------|-------------|------------------------|")
    for pid in ["P_H1", "P_H2", "P_H3", "P_H4"]:
        p = predictions[pid]
        # escape pipes in empirical text
        emp = p["empirical_value"].replace("|", "\\|")
        lines.append(f"| {pid}  | **{p['verdict']}** | {emp} |")
    lines.append("")

    lines.append("## Prediction statements")
    lines.append("")
    for pid in ["P_H1", "P_H2", "P_H3", "P_H4"]:
        p = predictions[pid]
        lines.append(f"### {pid} -- {p['verdict']}")
        lines.append("")
        lines.append(f"**Statement:** {p['statement']}")
        lines.append("")
        lines.append(f"**Empirical:** {p['empirical_value']}")
        lines.append("")
        lines.append(f"**Notes:** {p['notes']}")
        lines.append("")

    lines.append("## Crossover generation per scenario")
    lines.append("")
    lines.append("| Scenario | Dominant mode | Reps reaching >95% | Median crossover gen | Range |")
    lines.append("|----------|---------------|---------------------|-----------------------|--------|")
    for scen in ["A", "B", "C", "D"]:
        c = crossovers[scen]
        rng_str = f"{c['min_crossover_gen']}--{c['max_crossover_gen']}" if c['n_reps_reached_95pct'] > 0 else "--"
        lines.append(f"| {scen} | {c['dominant_mode']} | "
                     f"{c['n_reps_reached_95pct']}/{c['n_reps_total']} | "
                     f"{c['median_crossover_gen']} | {rng_str} |")
    lines.append("")

    lines.append("## Per-scenario final frequency (mean across 10 replicates)")
    lines.append("")
    lines.append("| Scenario | Mode | Mean final freq | SD final freq | Mean final fitness |")
    lines.append("|----------|------|------------------|----------------|---------------------|")
    for row in scenario_summary_table:
        lines.append(f"| {row['scenario']} | {row['mode_label']} | "
                     f"{row['mean_final_freq']:.4f} | {row['sd_final_freq']:.4f} | "
                     f"{row['mean_final_fitness']:.4f} |")
    lines.append("")

    lines.append("## Reproducibility")
    lines.append("")
    lines.append("- Module seed: `np.random.seed(42)` and `np.random.default_rng(42)`")
    lines.append("- Per-replicate seeds: `42 + r` for `r in {0..9}`")
    lines.append("- Mode 6 fixed template per replicate: derived from the per-replicate seed")
    lines.append(f"- K = {K}, L_TARGET = {L_TARGET}, N_GEN = {N_GEN}, BETA = {BETA}")
    lines.append(f"- MU_MODE1 = {MU_MODE1}, MU_MODE5 = {MU_MODE5}, N_MODULES_5 = {N_MODULES_5}")
    lines.append(f"- EPS_NOISE_5 = {EPS_NOISE_5}, EPS_NOISE_6 = {EPS_NOISE_6}")
    lines.append(f"- 10 replicates per scenario, 4 scenarios, {N_GEN} generations each")
    lines.append("")
    lines.append("## Mode 6 implementation note")
    lines.append("")
    lines.append("Implementation A (binary copyability) is the primary test. Mode 6 has a single fixed "
                 "template per replicate (the 2D-surface pattern), and every Mode 6 agent draws its "
                 "phenotype each generation as `template + i.i.d. noise`. There is no copy mechanism, so "
                 "selection within Mode 6 cannot move the population mean phenotype away from the fixed "
                 "template. Implementations B (per-lineage fixed template) and C (heritability < 1) "
                 "from dispatch3.md §Step 4 are deferred to a sensitivity follow-up.")
    lines.append("")
    Path(path).write_text("\n".join(lines))


def buildScenarioSummaryTable(summary_rows):
    out = []
    by_scen_mode = {}
    for r in summary_rows:
        by_scen_mode.setdefault((r["scenario"], r["mode_label"]), []).append(r)
    for (scen, mode), rs in sorted(by_scen_mode.items()):
        ffs = [r["final_frequency"] for r in rs]
        ffts = [r["final_mean_fitness"] for r in rs if not (r["final_mean_fitness"] is None
                                                            or np.isnan(r["final_mean_fitness"]))]
        out.append({
            "scenario": scen,
            "mode_label": mode,
            "mean_final_freq": float(np.mean(ffs)),
            "sd_final_freq": float(np.std(ffs)),
            "mean_final_fitness": float(np.mean(ffts)) if ffts else float("nan"),
        })
    return out


# ============================================================================
# main
# ============================================================================
def main():
    np.random.seed(42)  # honor the project convention even though we don't use module-level draws

    progress_path = RESULTS_DIR / "test_h_progress.txt"
    progress_f = open(progress_path, "a", buffering=1)
    progress_f.write(f"\n# === test_h run started at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")

    t_start = time.time()

    # ---- Scenario A: Mode 1 (200) vs Mode 6 (200) ----
    progress_f.write("\n--- scenario A: Mode 1 (200) vs Mode 6 (200) ---\n")
    rows_a = runScenarioReplicates("A", makeMode1, makeMode6_factory, 200, 200,
                                    uses_mode6=True, progress_writer=progress_f)
    writeScenarioCsv(rows_a, "A", RESULTS_DIR / "test_h_scenario_A_v1.csv")
    plotScenarioFreqVsGen(rows_a, "A", FIGURES_DIR / "test_h_scenario_A_freq_vs_gen.png")

    # ---- Scenario B: Mode 1 (40) vs Mode 6 (360) ----
    progress_f.write("\n--- scenario B: Mode 1 (40) vs Mode 6 (360) ---\n")
    rows_b = runScenarioReplicates("B", makeMode1, makeMode6_factory, 40, 360,
                                    uses_mode6=True, progress_writer=progress_f)
    writeScenarioCsv(rows_b, "B", RESULTS_DIR / "test_h_scenario_B_v1.csv")

    # ---- Scenario C: Mode 5 (200) vs Mode 6 (200) ----
    progress_f.write("\n--- scenario C: Mode 5 (200) vs Mode 6 (200) ---\n")
    rows_c = runScenarioReplicates("C", makeMode5, makeMode6_factory, 200, 200,
                                    uses_mode6=True, progress_writer=progress_f)
    writeScenarioCsv(rows_c, "C", RESULTS_DIR / "test_h_scenario_C_v1.csv")

    # ---- Scenario D: Mode 1 (200) vs Mode 5 (200) ----
    progress_f.write("\n--- scenario D: Mode 1 (200) vs Mode 5 (200) ---\n")
    rows_d = runScenarioReplicates("D", makeMode1, makeMode5, 200, 200,
                                    uses_mode6=False, progress_writer=progress_f)
    writeScenarioCsv(rows_d, "D", RESULTS_DIR / "test_h_scenario_D_v1.csv")

    all_rows = rows_a + rows_b + rows_c + rows_d

    # ---- summary ----
    summary_rows = []
    summary_rows += summarizeScenario(all_rows, "A", "Mode1", "Mode6_implA")
    summary_rows += summarizeScenario(all_rows, "B", "Mode1", "Mode6_implA")
    summary_rows += summarizeScenario(all_rows, "C", f"Mode5_N{N_MODULES_5}", "Mode6_implA")
    summary_rows += summarizeScenario(all_rows, "D", "Mode1", f"Mode5_N{N_MODULES_5}")
    writeSummaryCsv(summary_rows, RESULTS_DIR / "test_h_summary_v1.csv")

    # ---- multi-scenario figure ----
    plotAllScenarios(all_rows, FIGURES_DIR / "test_h_all_scenarios.png")

    # ---- predictions ----
    predictions = evaluatePredictions(summary_rows)
    crossovers = crossoverPerScenario(summary_rows)
    summary_table = buildScenarioSummaryTable(summary_rows)
    writePredictionsMd(predictions, crossovers, summary_table,
                       RESULTS_DIR / "test_h_predictions_v1.md")

    # ---- print report ----
    dt_total = time.time() - t_start
    print()
    print("=" * 90)
    print(f"Test H finished in {dt_total/60:.1f} min")
    print("=" * 90)
    print()
    print("Pre-registered prediction results:")
    for pid in ["P_H1", "P_H2", "P_H3", "P_H4"]:
        p = predictions[pid]
        print(f"  {pid}: {p['verdict']}")
        print(f"    {p['empirical_value']}")
    print()
    print("Crossover generations per scenario (median across 10 reps):")
    for scen in ["A", "B", "C", "D"]:
        c = crossovers[scen]
        print(f"  {scen}: dominant={c['dominant_mode']}  "
              f"reached_95pct={c['n_reps_reached_95pct']}/{c['n_reps_total']}  "
              f"median_gen={c['median_crossover_gen']}  "
              f"range={c['min_crossover_gen']}-{c['max_crossover_gen']}")
    print()
    print(f"  predictions written: {RESULTS_DIR / 'test_h_predictions_v1.md'}")
    print("=" * 90)

    progress_f.write(f"\n# === test_h finished at {time.strftime('%Y-%m-%d %H:%M:%S')} (dt={dt_total:.1f}s) ===\n")
    progress_f.close()


if __name__ == "__main__":
    main()
