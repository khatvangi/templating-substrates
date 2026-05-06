"""
Test F -- family-level Mode 3 prediction sweep across the Drt3b family alignment.

Generalizes Test E v2 (one biological anchor: WT + E26Q + AbiK) to ~10 natural
clade variants present in the 1,232-sequence Deng et al. 2026 alignment.

Framework claim being tested
----------------------------
The Templating Substrates Framework partitions the Drt3b active-site residues:

  - UNIVERSAL GATES   R253 (100% conserved), G248 (99.9% conserved)
                      enforce the cycle architecture (Mode 3, N=2).
  - PRIMARY GATES     E26 (90%), Y289 (94%), Y650 (91%)
                      enforce per-state nucleotide selectivity.
  - SECONDARY RESIDUES R168, R408, Y170, T335, T338 (highly variable)
                      modulate fidelity epsilon WITHOUT altering N.

Mapped to the test_e_v2 channel parameterization, the framework predicts:

  - I_struct     ~ 1 - H_2(epsilon) bits per site, saturating at 1 bit
                  for clades with intact primary gates.
  - marginal dG  determined entirely by E26 identity (state_A selectivity):
                    E26 = E -> ~0.003   (WT-like)
                    E26 = D -> ~0.10    (intermediate, carboxylate retained
                                         but geometry shifted; analytically
                                         set to half-broken state_A)
                    E26 = Q -> ~0.10    (broken state_A -- Test E v2 result;
                                         dG misincorp 0.20 at state_A, so
                                         marginal G = 0.5 * 0.20 = 0.10)
                    E26 = other -> ~0.10 (treat as E26Q-like)
  - period-2 peak ~ (1-eps)^2 + eps^2; intact gates ~0.98, degraded ~0.83.
  - separation ratio ~234-1031x intact, ~100x for E26Q-degraded.

The strongest claim: secondary-residue-only substitutions DO NOT shift
predicted observables outside the Mode 3 N=2 envelope. Any clade observed
outside the envelope (e.g. marginal G > 0.30 or I_struct < 0.5 bits at a
secondary-only-variant clade) falsifies the framework's residue partitioning.

Parameterization rule (this rule IS the prediction)
---------------------------------------------------
For each clade signature (E26, R168, Y170, G248, R253, Y289, T335, T338, R408, Y650):

  state_A channel parameters (controls A,G selection at intended-A cycle phase):
    if E26 == "E":  P(A)=0.99,  P(G)=0.0033, P(C)=0.0033, P(T)=0.0034    eps_A=0.01
    elif E26 == "Q": P(A)=0.80, P(G)=0.20,  P(C)=0.0001, P(T)=0.0001    eps_A=0.20
                     (matches Test E v2 calibration to Deng et al. 2026)
    elif E26 == "D": P(A)=0.85, P(G)=0.10,  P(C)=0.025, P(T)=0.025      eps_A=0.15
                     (carboxylate retained but geometry shifted; halfway
                      between E and Q on the dG-misincorp axis)
    else:           same as Q   (treat all non-EDQ as broken state_A)

  state_C channel parameters (controls C,T selection at intended-C cycle phase):
    if Y289 == "Y" and G248 == "G":
        P(C)=0.99, P(A)=0.0033, P(G)=0.0033, P(T)=0.0034   eps_C=0.01
    elif G248 != "G":
        # broken steric exclusion of dG
        P(C)=0.85, P(G)=0.10, P(A)=0.025, P(T)=0.025       eps_C=0.15
    elif Y289 != "Y":
        # broken pyrimidine recognition
        P(C)=0.80, P(T)=0.15, P(A)=0.025, P(G)=0.025       eps_C=0.20
    else: unreachable

  cycle architecture:
    if R253 != "R":  cycle disrupted; predict I_struct < 0.5 bits.
                     model: collapse to single mixed channel (no phase coupling)
    else:            two-state cyclic Markov chain (Mode 3, N=2)

  Y650 (C-terminal protein-priming Tyr): does not enter per-state channel,
        but if absent we mark the clade as initiation-impaired (separate flag,
        does not change predicted observables for elongation).

  All five SECONDARY residues (R168, R408, Y170, T335, T338): no parameter change.
        The framework prediction is that their identity does NOT shift any
        predicted observable. CSV column 'envelope_classification' tracks this.

PASS / falsifiability
---------------------
PASS = secondary-residue-only clades all stay within the WT-like envelope:
   marginal G < 0.05, period-2 peak > 0.93, I_struct > 0.95 bits at L=64.
   AND primary-gate-degraded clades (E26 != E) shift in the predicted
   directions (marginal G > 0.05, period-2 peak < 0.93).
"""

import csv
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"
ALIGNMENT_PATH = PROJECT_ROOT / "science.aed1656_data_s1_to_s8 deng" / "science.aed1656_data_s3.fa"

REF_ACCESSION = "WP_126681219.1"

# (residue_position_1based, expected_wt_aa, role, class)
# class in {"universal", "primary", "secondary"}
TARGETS = [
    (26,  "E", "gatekeeper for dA selection",                "primary"),
    (168, "R", "pyrimidine-specific contact at dC15",        "secondary"),
    (170, "Y", "non-specific contact at dA14",               "secondary"),
    (248, "G", "steric exclusion of dG",                     "universal"),
    (253, "R", "cation-pi for dA / H-bond with dC17",        "universal"),
    (289, "Y", "pyrimidine-specific contact at dC15",        "primary"),
    (335, "T", "purine-specific contact at dA16",            "secondary"),
    (338, "T", "purine-specific contact at dA16",            "secondary"),
    (408, "R", "base-specific contact at dC13",              "secondary"),
    (650, "Y", "C-terminal protein-priming Tyr",             "primary"),
]

POSITIONS = [t[0] for t in TARGETS]
WT_AAS = [t[1] for t in TARGETS]
RESIDUE_CLASS = {t[0]: t[3] for t in TARGETS}
WT_TUPLE = tuple(WT_AAS)

PRIMARY_POSITIONS = [t[0] for t in TARGETS if t[3] == "primary"]
SECONDARY_POSITIONS = [t[0] for t in TARGETS if t[3] == "secondary"]
UNIVERSAL_POSITIONS = [t[0] for t in TARGETS if t[3] == "universal"]

# ---------------------------------------------------------------------------
# step 1: parse alignment and build per-sequence residue tuples
# ---------------------------------------------------------------------------

def parseAlignedFasta(path):
    """walk the file once, return list of (header, aligned_seq) tuples."""
    records = []
    header = None
    chunks = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            records.append((header, "".join(chunks)))
    return records


def buildPositionToColumnMap(aligned_seq, target_positions):
    """walk aligned_seq, count non-gap chars, return {pos -> column_index}."""
    wanted = set(target_positions)
    mapping = {}
    residue_count = 0
    for col_idx, aa in enumerate(aligned_seq):
        if aa == "-":
            continue
        residue_count += 1
        if residue_count in wanted:
            mapping[residue_count] = col_idx
            if len(mapping) == len(wanted):
                break
    return mapping


def extractResidueTuple(aligned_seq, columns):
    """pull aa at each column; gap or out-of-range -> '-'."""
    out = []
    for c in columns:
        if c >= len(aligned_seq):
            out.append("-")
        else:
            out.append(aligned_seq[c])
    return tuple(out)


def findClades(records, columns, min_members=5):
    """group sequences by 10-residue signature, drop tuples with n < min_members.

    returns list of dicts sorted by descending member count.
    """
    bag = defaultdict(list)
    for hdr, seq in records:
        sig = extractResidueTuple(seq, columns)
        # drop signatures with any gap/X (cannot make a confident prediction)
        if "-" in sig or "X" in sig:
            continue
        bag[sig].append(hdr)
    clades = []
    for sig, members in bag.items():
        if len(members) >= min_members:
            clades.append({"signature": sig, "members": members, "n": len(members)})
    clades.sort(key=lambda c: -c["n"])
    return clades


# ---------------------------------------------------------------------------
# step 2: parameterization rule -- maps clade signature to channel matrices
# ---------------------------------------------------------------------------
# rule for each residue is documented in the module docstring.

WT_STATE_A = np.array([0.99, 0.0033, 0.0033, 0.0034])
WT_STATE_C = np.array([0.0033, 0.99, 0.0033, 0.0034])

E26Q_STATE_A = np.array([0.80, 0.0001, 0.20, 0.0001])  # P(A,C,G,T) -- Deng et al. 2026
E26D_STATE_A = np.array([0.85, 0.025, 0.10, 0.025])    # intermediate, half-broken

G248_BROKEN_STATE_C = np.array([0.025, 0.85, 0.10, 0.025])   # broken steric exclusion of dG
Y289_BROKEN_STATE_C = np.array([0.025, 0.80, 0.025, 0.15])   # broken pyrimidine recognition

# alphabet index mapping for clarity
A, C, G, T = 0, 1, 2, 3


def findClassesForClade(signature):
    """compute (state_A_p, state_C_p, cycle_intact, init_intact) from signature."""
    sig = dict(zip(POSITIONS, signature))
    e26  = sig[26]
    g248 = sig[248]
    r253 = sig[253]
    y289 = sig[289]
    y650 = sig[650]

    # state_A determined by E26 identity
    if e26 == "E":
        state_A = WT_STATE_A
    elif e26 == "Q":
        state_A = E26Q_STATE_A
    elif e26 == "D":
        state_A = E26D_STATE_A
    else:
        # treat all other substitutions as E26Q-like (broken state_A)
        state_A = E26Q_STATE_A

    # state_C determined by G248 then Y289
    if g248 == "G" and y289 == "Y":
        state_C = WT_STATE_C
    elif g248 != "G":
        state_C = G248_BROKEN_STATE_C
    else:
        state_C = Y289_BROKEN_STATE_C

    cycle_intact = (r253 == "R")
    init_intact  = (y650 == "Y" or y650 == "F")  # Y/F preserves OH/aromatic
    return state_A, state_C, cycle_intact, init_intact


def envelopeClassification(signature):
    """tag clade based on which residues differ from WT."""
    sig = dict(zip(POSITIONS, signature))
    primary_diff = [p for p in PRIMARY_POSITIONS
                    if sig[p] != dict(zip(POSITIONS, WT_TUPLE))[p]]
    universal_diff = [p for p in UNIVERSAL_POSITIONS
                      if sig[p] != dict(zip(POSITIONS, WT_TUPLE))[p]]
    secondary_diff = [p for p in SECONDARY_POSITIONS
                      if sig[p] != dict(zip(POSITIONS, WT_TUPLE))[p]]

    if universal_diff:
        return "out-of-envelope (universal-gate disruption)"
    if primary_diff:
        return f"shifted (primary-gate substitution at {primary_diff})"
    if secondary_diff:
        return "in-envelope (secondary-only variation)"
    return "in-envelope (WT)"


# ---------------------------------------------------------------------------
# step 3: simulator and estimator -- duplicated verbatim from test_e_v2
# (per repo CLAUDE.md rule 4: tests are self-contained on purpose)
# ---------------------------------------------------------------------------

def _draw_from(p, size, rng):
    """vectorized categorical draw from probability vector p."""
    cum = np.cumsum(p)
    u = rng.random(size=size)
    return np.searchsorted(cum, u).astype(np.int16)


def simulate_two_state_cycle(L, n_samples, state_A_p, state_C_p, rng):
    """generic two-state cyclic active-site simulator."""
    X = rng.integers(0, 2, size=n_samples, dtype=np.int16)  # phase
    pos = np.arange(L, dtype=np.int64)[None, :]
    intended = (X[:, None].astype(np.int64) + pos) % 2  # 0 = state_A, 1 = state_C
    yA = _draw_from(state_A_p, (n_samples, L), rng)
    yC = _draw_from(state_C_p, (n_samples, L), rng)
    Y = np.where(intended == 0, yA, yC).astype(np.int16)
    return X, Y


def simulate_collapsed_no_cycle(L, n_samples, state_A_p, state_C_p, rng):
    """fallback when R253 is mutated: no cycle, every position is the
    average channel; X exists but has no causal effect (predicts low I_struct)."""
    avg = 0.5 * (state_A_p + state_C_p)
    X = rng.integers(0, 2, size=n_samples, dtype=np.int16)
    Y = _draw_from(avg, (n_samples, L), rng).astype(np.int16)
    return X, Y


EPS_CLIP = 1e-12


def _entropy_from_counts(counts):
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts.astype(np.float64) / total
    nz = p[p > 0]
    return float(-(nz * np.log2(nz)).sum())


def estimate_joint_mi(X, Y, x_alphabet=2, y_alphabet=4, block_length=None):
    """plug-in joint MI: I(X; Y[0:k]) = H(Y[0:k]) - H(Y[0:k] | X)."""
    n_samples, L = Y.shape
    if block_length is None:
        block_length = min(L, 8)
    k = min(block_length, L)
    if k < 1:
        return 0.0, 0, 0
    while y_alphabet ** k * x_alphabet > n_samples / 10 and k > 1:
        k -= 1
    powers = y_alphabet ** np.arange(k, dtype=np.int64)
    y_codes = (Y[:, :k].astype(np.int64) * powers[None, :]).sum(axis=1)
    counts_y = np.bincount(y_codes, minlength=y_alphabet ** k)
    H_y = _entropy_from_counts(counts_y)
    H_y_cond = 0.0
    counts_x = np.bincount(X.astype(np.int64), minlength=x_alphabet)
    for x in range(x_alphabet):
        nx = counts_x[x]
        if nx == 0:
            continue
        mask = X == x
        cyx = np.bincount(y_codes[mask], minlength=y_alphabet ** k)
        H_y_cond += (nx / n_samples) * _entropy_from_counts(cyx)
    return H_y - H_y_cond, k, H_y


def estimate_periodicity(Y, max_lag):
    """symbol-equality autocorrelation."""
    _, L = Y.shape
    autocorr = np.zeros(max_lag, dtype=np.float64)
    for lag in range(1, max_lag + 1):
        if lag >= L:
            autocorr[lag - 1] = 0.0
            continue
        eq = (Y[:, : L - lag] == Y[:, lag:])
        autocorr[lag - 1] = float(eq.mean())
    return autocorr


def simulate_bulk_matched(L, n_samples, target_marginal, x_alphabet, rng):
    X_bulk = rng.integers(0, x_alphabet, size=n_samples, dtype=np.int16)
    Y_bulk = _draw_from(target_marginal, (n_samples, L), rng)
    return X_bulk, Y_bulk


def empirical_marginal(Y, alphabet_size=4):
    counts = np.bincount(Y.ravel(), minlength=alphabet_size)
    return counts / counts.sum()


# ---------------------------------------------------------------------------
# step 4: per-clade simulation driver
# ---------------------------------------------------------------------------

def runSimsForClade(state_A, state_C, cycle_intact, L_values, n_reps, n_samples, rng):
    """run mode-3 sim at each L for n_reps reps; return per-L stats."""
    sim_fn = simulate_two_state_cycle if cycle_intact else simulate_collapsed_no_cycle
    out = []
    for L in L_values:
        I_emp_reps, I_bulk_reps, ratio_reps = [], [], []
        peak_lag_reps, peak_val_reps = [], []
        margA, margC, margG, margT = [], [], [], []
        for _ in range(n_reps):
            X, Y = sim_fn(L, n_samples, state_A, state_C, rng)
            target_marg = empirical_marginal(Y)
            Xb, Yb = simulate_bulk_matched(L, n_samples, target_marg, 2, rng)
            I_emp, _, _ = estimate_joint_mi(X, Y)
            I_bulk, _, _ = estimate_joint_mi(Xb, Yb)
            ratio = I_emp / max(I_bulk, EPS_CLIP)
            max_lag = min(8, L - 1) if L > 1 else 1
            if max_lag >= 1:
                ac = estimate_periodicity(Y, max_lag=max_lag)
                # specifically look at lag-2 (predicted period of Mode 3 N=2)
                # but also report the absolute peak
                peak_idx = int(np.argmax(ac))
                peak_lag_reps.append(peak_idx + 1)
                peak_val_reps.append(float(ac[peak_idx]))
            else:
                peak_lag_reps.append(0)
                peak_val_reps.append(0.0)
            I_emp_reps.append(I_emp)
            I_bulk_reps.append(I_bulk)
            ratio_reps.append(ratio)
            margA.append(target_marg[A]); margC.append(target_marg[C])
            margG.append(target_marg[G]); margT.append(target_marg[T])
        # explicit lag-2 peak (the framework's predicted period)
        # re-derive once on a single batch for the lag-2 number
        X1, Y1 = sim_fn(L, n_samples, state_A, state_C, rng)
        if L > 2:
            ac1 = estimate_periodicity(Y1, max_lag=min(8, L - 1))
            lag2_val = float(ac1[1])
        else:
            lag2_val = 0.0
        out.append({
            "L": L,
            "n_reps": n_reps,
            "n_samples_per_rep": n_samples,
            "I_struct_mean": float(np.mean(I_emp_reps)),
            "I_struct_std": float(np.std(I_emp_reps)),
            "I_bulk_mean": float(np.mean(I_bulk_reps)),
            "separation_ratio_mean": float(np.mean(ratio_reps)),
            "peak_lag_mode": int(Counter(peak_lag_reps).most_common(1)[0][0]),
            "peak_val_mean": float(np.mean(peak_val_reps)),
            "lag2_peak_value": lag2_val,
            "marg_A_mean": float(np.mean(margA)),
            "marg_C_mean": float(np.mean(margC)),
            "marg_G_mean": float(np.mean(margG)),
            "marg_T_mean": float(np.mean(margT)),
        })
    return out


# ---------------------------------------------------------------------------
# step 5: build prediction table from per-clade sims
# ---------------------------------------------------------------------------

def summarizeClade(clade, sim_rows, ref_L=64):
    """pull headline numbers from L=ref_L (or nearest)."""
    chosen = min(sim_rows, key=lambda r: abs(r["L"] - ref_L))
    return {
        "predicted_marginal_G_atL": chosen["marg_G_mean"],
        "predicted_I_struct_atL": chosen["I_struct_mean"],
        "predicted_period2_peak_atL": chosen["lag2_peak_value"],
        "predicted_separation_ratio_atL": chosen["separation_ratio_mean"],
        "summary_L": chosen["L"],
    }


# ---------------------------------------------------------------------------
# main driver
# ---------------------------------------------------------------------------

def shortenHeader(hdr):
    # extract accession-like token
    tok = hdr.split()[0].lstrip(">")
    return tok[:40]


def main():
    np.random.seed(42)
    rng = np.random.default_rng(42)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    progress_path = RESULTS_DIR / "test_f_progress.txt"
    pf = open(progress_path, "w", buffering=1)

    pf.write(f"# === test_f run started at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
    print(f"[info] parsing alignment {ALIGNMENT_PATH}")
    records = parseAlignedFasta(str(ALIGNMENT_PATH))
    n_seqs = len(records)
    print(f"[info] parsed {n_seqs} sequences")
    pf.write(f"[info] parsed {n_seqs} sequences\n")

    # find ref row, build pos -> col map
    ref_idx = next((i for i, (h, _) in enumerate(records) if REF_ACCESSION in h), None)
    if ref_idx is None:
        raise SystemExit(f"reference {REF_ACCESSION} not found in alignment")
    _, ref_aligned = records[ref_idx]
    pos_to_col = buildPositionToColumnMap(ref_aligned, POSITIONS)
    columns = [pos_to_col[p] for p in POSITIONS]
    print(f"[info] mapped positions {POSITIONS} -> columns {columns}")
    pf.write(f"[info] columns = {columns}\n")

    # find clades (n >= 5)
    clades = findClades(records, columns, min_members=5)
    print(f"[info] found {len(clades)} clades with n >= 5")
    pf.write(f"[info] found {len(clades)} clades with n >= 5\n")

    # assign IDs and labels
    for i, c in enumerate(clades, start=1):
        c["clade_id"] = f"C{i:02d}"
        if c["signature"] == WT_TUPLE:
            c["clade_id"] += "_WT"
        # tag by primary-gate substitution for human readability
        sig = dict(zip(POSITIONS, c["signature"]))
        prim_diffs = [(p, sig[p]) for p in PRIMARY_POSITIONS if sig[p] != WT_TUPLE[POSITIONS.index(p)]]
        univ_diffs = [(p, sig[p]) for p in UNIVERSAL_POSITIONS if sig[p] != WT_TUPLE[POSITIONS.index(p)]]
        if univ_diffs:
            c["clade_id"] += "_UNIV-" + ",".join(f"{wt}{p}{aa}" for (p, aa), wt in
                                                   zip(univ_diffs, [WT_TUPLE[POSITIONS.index(p)] for p, _ in univ_diffs]))
        elif prim_diffs:
            c["clade_id"] += "_PRIM-" + ",".join(f"{wt}{p}{aa}" for (p, aa), wt in
                                                   zip(prim_diffs, [WT_TUPLE[POSITIONS.index(p)] for p, _ in prim_diffs]))

    # report top clades
    print("\n[info] top 10 clades by member count:")
    for c in clades[:10]:
        print(f"  {c['clade_id']:<40s}  n={c['n']:>5d}  sig={''.join(c['signature'])}")
        pf.write(f"clade {c['clade_id']}  n={c['n']}  sig={''.join(c['signature'])}\n")

    # parameters
    L_values   = [4, 8, 16, 32, 64, 128, 256, 500]
    n_reps     = 30      # gives reasonable error bars
    n_samples  = 5000    # per rep; 30 * 5000 = 150,000 effective draws per L

    # run sims per clade
    per_clade_rows = []
    prediction_rows = []
    print("\n[info] running mode-3 sims per clade...")
    print(f"{'clade_id':<40s} {'n':>5s} {'classification':<55s} {'eta(s)':>7s}")
    print("-" * 115)

    t_start_all = time.time()
    for c in clades:
        t0 = time.time()
        state_A, state_C, cycle_intact, init_intact = findClassesForClade(c["signature"])
        sim_rows = runSimsForClade(state_A, state_C, cycle_intact,
                                    L_values, n_reps, n_samples, rng)
        for sr in sim_rows:
            row = {"clade_id": c["clade_id"],
                   "signature": "".join(c["signature"]),
                   "n_members": c["n"],
                   "cycle_intact": int(cycle_intact),
                   "init_intact": int(init_intact)}
            row.update(sr)
            per_clade_rows.append(row)

        summary = summarizeClade(c, sim_rows, ref_L=64)
        env_class = envelopeClassification(c["signature"])
        prediction_rows.append({
            "clade_id": c["clade_id"],
            "signature": "".join(c["signature"]),
            "n_members": c["n"],
            "representative": shortenHeader(c["members"][0]),
            "cycle_intact": int(cycle_intact),
            "init_intact": int(init_intact),
            "predicted_marginal_G": summary["predicted_marginal_G_atL"],
            "predicted_I_struct": summary["predicted_I_struct_atL"],
            "predicted_period2_peak": summary["predicted_period2_peak_atL"],
            "predicted_separation_ratio": summary["predicted_separation_ratio_atL"],
            "summary_L": summary["summary_L"],
            "envelope_classification": env_class,
        })
        elapsed = time.time() - t0
        print(f"{c['clade_id']:<40s} {c['n']:>5d} {env_class:<55s} {elapsed:>7.1f}")
        pf.write(f"clade {c['clade_id']}: predicted G={summary['predicted_marginal_G_atL']:.4f}, "
                 f"I={summary['predicted_I_struct_atL']:.4f}, "
                 f"per2={summary['predicted_period2_peak_atL']:.4f}, "
                 f"sep={summary['predicted_separation_ratio_atL']:.1f}x, "
                 f"elapsed={elapsed:.1f}s\n")
    total_elapsed = time.time() - t_start_all
    print(f"\n[info] all sims finished in {total_elapsed:.1f}s")
    pf.write(f"[info] total sim time = {total_elapsed:.1f}s\n")
    pf.close()

    # write CSVs
    save_predictions_csv(prediction_rows, RESULTS_DIR / "test_f_family_predictions_v1.csv")
    save_per_clade_csv(per_clade_rows, RESULTS_DIR / "test_f_per_clade_simulations_v1.csv")

    # falsifiability statement
    write_falsifiability(prediction_rows, RESULTS_DIR / "test_f_falsifiability_v1.md")

    # figure
    plot_clade_predictions(prediction_rows, FIGURES_DIR / "test_f_clade_predictions.png")

    # readme
    write_readme(RESULTS_DIR / "test_f_README.md", prediction_rows, len(clades),
                 total_elapsed, n_seqs)

    # PASS/FAIL evaluation
    verdict = evaluatePassFail(prediction_rows)
    print("\n" + "=" * 95)
    print(f"OVERALL: {verdict['overall']}")
    print("=" * 95)
    for line in verdict["lines"]:
        print(line)


def save_predictions_csv(rows, path):
    fields = ["clade_id", "signature", "n_members", "representative",
              "cycle_intact", "init_intact",
              "predicted_marginal_G", "predicted_I_struct",
              "predicted_period2_peak", "predicted_separation_ratio",
              "summary_L", "envelope_classification"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def save_per_clade_csv(rows, path):
    fields = ["clade_id", "signature", "n_members",
              "cycle_intact", "init_intact",
              "L", "n_reps", "n_samples_per_rep",
              "I_struct_mean", "I_struct_std",
              "I_bulk_mean", "separation_ratio_mean",
              "peak_lag_mode", "peak_val_mean", "lag2_peak_value",
              "marg_A_mean", "marg_C_mean", "marg_G_mean", "marg_T_mean"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def write_falsifiability(rows, path):
    """compute X (max marginal G), Y (min period-2 peak), Z (min I_struct) over
    secondary-only clades; report explicit numbers."""
    sec_only = [r for r in rows if r["envelope_classification"].startswith("in-envelope")]
    primary  = [r for r in rows if r["envelope_classification"].startswith("shifted")]
    universal = [r for r in rows if r["envelope_classification"].startswith("out-of-envelope")]

    if sec_only:
        X = max(r["predicted_marginal_G"] for r in sec_only)
        Y = min(r["predicted_period2_peak"] for r in sec_only)
        Z = min(r["predicted_I_struct"] for r in sec_only)
        K = len(sec_only)
    else:
        X, Y, Z, K = float("nan"), float("nan"), float("nan"), 0

    if primary:
        Gp = [r["predicted_marginal_G"] for r in primary]
        Pp = [r["predicted_period2_peak"] for r in primary]
        Kprime = len(primary)
        Xp_min = min(Gp); Xp_max = max(Gp)
        Yp_min = min(Pp); Yp_max = max(Pp)
    else:
        Kprime = 0
        Xp_min = Xp_max = Yp_min = Yp_max = float("nan")

    text = f"""# Test F -- falsifiability statement (v1)

## Predictions

Test F predicts that for the **K = {K}** largest secondary-residue clades retaining
intact universal gates (R253, G248) and primary gates (E26, Y289, Y650),
cDIP-seq on the natural family representatives would yield:

- marginal G fraction <= **X = {X:.4f}** at L=64
- period-2 autocorrelation peak >= **Y = {Y:.4f}** at L=64
- $I_\\text{{struct}}$ >= **Z = {Z:.4f}** bits at L=64

(Numbers above are the worst-case observable across all in-envelope clades
in the simulation; the framework predicts every member clade meets each
threshold.)

Any clade observed outside this envelope -- e.g. a secondary-residue-only
representative for which marginal G > {X:.3f} or $I_\\text{{struct}}$ < {Z:.3f}
bits -- **falsifies** the framework's claim that secondary residues
(R168, R408, Y170, T335, T338) do not affect mode classification.

Conversely, the **K' = {Kprime}** clades with naturally varying primary gates
predict shifted observables in the directions documented in the
parameterization rule (see code/test_f_family_sweep.py module docstring):

- marginal G in the simulated range [{Xp_min:.4f}, {Xp_max:.4f}]
- period-2 peak in the simulated range [{Yp_min:.4f}, {Yp_max:.4f}]

Confirming the predicted **direction** for these clades (G shifts up,
period-2 peak shifts down for E26->Q/D substitutions) strengthens the
framework. Observing **perpendicular** shifts -- secondary-residue clades
shifting while primary-gate clades do not -- falsifies it.

## Universal-gate disrupted clades

{len(universal)} clades (out of {len(rows)} total) carry universal-gate
substitutions (R253 != R or G248 != G). Framework predicts these are
out-of-Mode-3 entirely: $I_\\text{{struct}}$ < 0.5 bits, period-2 peak ~0.25
(no cycle structure). These are the strongest pre-registered out-of-envelope
predictions.

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
| marginal G (in-envel.)  | <= {X:.4f}                     |
| period-2 peak (in-env.) | >= {Y:.4f}                     |
| I_struct (in-env.)      | >= {Z:.4f} bits                |
| marginal G (primary)    | in [{Xp_min:.4f}, {Xp_max:.4f}] |
| period-2 peak (primary) | in [{Yp_min:.4f}, {Yp_max:.4f}] |
"""
    path.write_text(text)


def evaluatePassFail(rows):
    """check the secondary-only clades stay in envelope and per-state-relevant
    primary-gate clades (E26, Y289, G248) shift in the predicted direction.

    note: Y650 is a primary residue but the parameterization rule explicitly
    documents that Y650 controls initiation, not the per-state elongation
    channel. so a Y650-only-substituted clade is predicted to stay in envelope
    for the observables we measure (marginal G, period-2 peak, I_struct).
    """
    # envelope thresholds (matching Test E v2 WT performance, padded for sim noise):
    G_thresh = 0.05    # marginal G must stay below 0.05 for in-envelope
    P_thresh = 0.93    # period-2 peak must stay above 0.93 for in-envelope
    I_thresh = 0.95    # I_struct must stay above 0.95 bits for in-envelope

    # partition primary-gate clades by which residue is changed
    sec_only   = [r for r in rows if r["envelope_classification"].startswith("in-envelope")]
    primary    = [r for r in rows if r["envelope_classification"].startswith("shifted")]

    # subdivide primary by changed residue: per-state relevant vs Y650 (initiation only)
    per_state_primary = []  # E26 / Y289 / G248 changed -> framework predicts shift
    init_only_primary = []  # only Y650 changed -> framework predicts NO shift in elongation
    wt_idx = {p: i for i, p in enumerate(POSITIONS)}
    for r in primary:
        sig = r["signature"]
        prim_diffs = [p for p in PRIMARY_POSITIONS if sig[wt_idx[p]] != WT_TUPLE[wt_idx[p]]]
        if prim_diffs == [650]:
            init_only_primary.append(r)
        else:
            per_state_primary.append(r)

    sec_violations = []
    for r in sec_only:
        if (r["predicted_marginal_G"] > G_thresh or
            r["predicted_period2_peak"] < P_thresh or
            r["predicted_I_struct"] < I_thresh):
            sec_violations.append(r["clade_id"])

    init_only_violations = []
    for r in init_only_primary:
        # framework predicts these stay in envelope for elongation observables
        if (r["predicted_marginal_G"] > G_thresh or
            r["predicted_period2_peak"] < P_thresh or
            r["predicted_I_struct"] < I_thresh):
            init_only_violations.append(r["clade_id"])

    per_state_correct_direction = []
    for r in per_state_primary:
        # per-state-relevant primary-gate clades should show G shifted up OR period-2 down
        if r["predicted_marginal_G"] > G_thresh or r["predicted_period2_peak"] < P_thresh:
            per_state_correct_direction.append(r["clade_id"])

    pass_sec        = len(sec_violations) == 0 and len(sec_only) > 0
    pass_init_only  = len(init_only_violations) == 0
    pass_per_state  = (len(per_state_primary) == 0) or \
                      (len(per_state_correct_direction) == len(per_state_primary))

    overall = "PASS" if (pass_sec and pass_init_only and pass_per_state) else "FAIL"
    lines = [
        f"Secondary-only clades:           {len(sec_only)} total, {len(sec_violations)} violations",
    ]
    if sec_violations:
        lines.append(f"  violators: {sec_violations}")
    lines.append(f"Per-state primary clades (E26/Y289/G248):  "
                 f"{len(per_state_primary)} total, "
                 f"{len(per_state_correct_direction)} shifted in predicted direction")
    lines.append(f"Initiation-only primary clades (Y650F): "
                 f"{len(init_only_primary)} total, {len(init_only_violations)} violations "
                 f"(framework predicts these stay in elongation envelope)")
    return {"overall": overall, "lines": lines,
            "pass_sec": pass_sec, "pass_per_state": pass_per_state,
            "pass_init_only": pass_init_only}


def plot_clade_predictions(rows, path):
    """scatter of predicted marginal G vs predicted period-2 peak, color by class."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    color_map = {
        "in-envelope (WT)": "#2563eb",
        "in-envelope (secondary-only variation)": "#10b981",
        "shifted (primary-gate substitution at [26])": "#f59e0b",
        "out-of-envelope (universal-gate disruption)": "#ef4444",
    }
    def colorOf(env_str):
        for k, v in color_map.items():
            if env_str.startswith(k.split(" (")[0]):
                if "WT" in env_str: return color_map["in-envelope (WT)"]
                if "secondary" in env_str: return color_map["in-envelope (secondary-only variation)"]
                if "shifted" in env_str: return color_map["shifted (primary-gate substitution at [26])"]
                if "out-of-envelope" in env_str: return color_map["out-of-envelope (universal-gate disruption)"]
        return "#6b7280"

    ax = axes[0]
    for r in rows:
        c = colorOf(r["envelope_classification"])
        size = 30 + 5 * np.log2(max(r["n_members"], 2))
        ax.scatter(r["predicted_marginal_G"], r["predicted_period2_peak"],
                   s=size, color=c, edgecolor="black", linewidth=0.5, alpha=0.85)
        ax.annotate(r["clade_id"].split("_")[0],
                    (r["predicted_marginal_G"], r["predicted_period2_peak"]),
                    fontsize=7, alpha=0.7,
                    xytext=(4, 3), textcoords="offset points")
    ax.set_xlabel("predicted marginal G fraction (at L=64)")
    ax.set_ylabel("predicted period-2 autocorrelation peak (at L=64)")
    ax.set_title("Test F: clade-level predictions in (marg G, period-2) plane")
    ax.axvline(0.05, color="grey", linestyle="--", alpha=0.4, label="G threshold = 0.05")
    ax.axhline(0.93, color="grey", linestyle=":", alpha=0.4, label="period-2 threshold = 0.93")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

    ax = axes[1]
    for r in rows:
        c = colorOf(r["envelope_classification"])
        size = 30 + 5 * np.log2(max(r["n_members"], 2))
        ax.scatter(r["predicted_I_struct"], r["predicted_separation_ratio"],
                   s=size, color=c, edgecolor="black", linewidth=0.5, alpha=0.85)
    ax.set_xlabel(r"predicted $I_{\mathrm{struct}}$ (bits, at L=64)")
    ax.set_ylabel("predicted separation ratio (vs bulk-matched control)")
    ax.set_yscale("log")
    ax.set_title("Test F: I_struct vs separation ratio per clade")
    ax.axvline(0.95, color="grey", linestyle="--", alpha=0.4, label="I_struct threshold = 0.95")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # hand-build legend across panels
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=v,
                          markeredgecolor="black", label=k.split(" (")[0])
               for k, v in color_map.items()]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=9,
               bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


README_TEMPLATE = """# Test F -- family-level Mode 3 prediction sweep (v1)

## What this tests

Generalizes Test E v2 (one biological anchor: Drt3b WT + E26Q + AbiK) to ~K
natural clade variants present in the 1,232-sequence Deng et al. 2026 alignment
of Drt3b homologs. Each clade is defined by its 10-residue signature at the
catalytic positions (E26, R168, Y170, G248, R253, Y289, T335, T338, R408, Y650).

## Parameterization rule (the prediction)

The framework's prediction is encoded entirely in the parameterization rule
mapping clade signature -> per-state channel. See module docstring of
`code/test_f_family_sweep.py` for the full rule. Headlines:

  - E26 = E         -> state_A intact (P(G)=0.0033, eps_A=0.01)
  - E26 = D         -> state_A half-broken (P(G)=0.10, eps_A=0.15)
  - E26 = Q / other -> state_A broken    (P(G)=0.20, eps_A=0.20)
  - G248 != G       -> state_C broken (dG misincorp)
  - Y289 != Y       -> state_C broken (pyrimidine recognition lost)
  - R253 != R       -> cycle disrupted; predict I_struct < 0.5 bits
  - secondary residues (R168, R408, Y170, T335, T338) -> NO parameter change
    (the framework predicts these do not shift observables)

## Outputs

  - test_f_family_predictions_v1.csv      one row per clade, with the
                                          framework's predicted observables
  - test_f_per_clade_simulations_v1.csv   raw L-sweep simulation outputs
  - test_f_falsifiability_v1.md           pre-registered numerical thresholds
  - figures/test_f_clade_predictions.png  scatter of predictions in observable space

## How to falsify

For any clade row in test_f_family_predictions_v1.csv, express the named
representative and run cDIP-seq. If a secondary-residue-only clade shifts
its observables outside the in-envelope thresholds (G > 0.05, period-2 < 0.93,
I_struct < 0.95), the framework's residue partitioning is wrong.

## Run summary

  - {n_seqs} sequences in alignment
  - {n_clades} clades found with n >= 5 members
  - sim time: {elapsed:.1f}s
  - L sweep: 4, 8, 16, 32, 64, 128, 256, 500
  - 30 reps x 5000 samples per (clade, L)

## Pass criterion

PASS = secondary-only clades all stay within the WT-like envelope at L=64
AND primary-gate-degraded clades all shift in the predicted direction.

See test_f_falsifiability_v1.md for explicit numerical thresholds.
"""


def write_readme(path, rows, n_clades, elapsed, n_seqs):
    txt = README_TEMPLATE.format(n_seqs=n_seqs, n_clades=n_clades, elapsed=elapsed)
    path.write_text(txt)


if __name__ == "__main__":
    main()
