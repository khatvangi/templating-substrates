"""
conservation analysis of paper-named gate residues in the Drt3b family alignment

input: aligned FASTA from Deng et al. 2026 supplementary data S3
reference: EcDrt3b (WP_126681219.1, 650 aa, starts MSKKKEVRVNKKDFNRVL)

steps:
  1. parse aligned FASTA (stdlib only)
  2. locate EcDrt3b row by accession in header
  3. sanity-check gap-stripped length and N-terminal residues
  4. map paper residue numbers (1-based) to alignment column indices
     by walking the EcDrt3b row and counting non-gap characters
  5. tally amino acids at each target column across all sequences
  6. emit a clean stdout table + CSV with top-5 amino acids per position
"""

import csv
import os
from collections import Counter

# ---------------------------------------------------------------------------
# config
# ---------------------------------------------------------------------------

ALIGNMENT_PATH = (
    "/storage/kiran-stuff/templating_framework/"
    "science.aed1656_data_s1_to_s8 deng/science.aed1656_data_s3.fa"
)
CSV_OUT = (
    "/storage/kiran-stuff/templating_framework/"
    "results/drt3b_conservation_analysis.csv"
)

REF_ACCESSION = "WP_126681219.1"
REF_EXPECTED_LEN = 650
REF_EXPECTED_PREFIX = "MSKKKEVRVNKKDFNRVL"

# (residue_position_1based, expected_wt_aa, role)
TARGETS = [
    (26,  "E", "gatekeeper for dA selection"),
    (168, "R", "pyrimidine-specific contact at dC15"),
    (170, "Y", "non-specific contact at dA14"),
    (248, "G", "steric exclusion of dG"),
    (253, "R", "cation-pi for dA / H-bond with dC17"),
    (289, "Y", "pyrimidine-specific contact at dC15"),
    (335, "T", "purine-specific contact at dA16"),
    (338, "T", "purine-specific contact at dA16"),
    (408, "R", "base-specific contact at dC13"),
    (650, "Y", "C-terminal protein-priming Tyr"),
]


# ---------------------------------------------------------------------------
# fasta parsing
# ---------------------------------------------------------------------------

def parseAlignedFasta(path):
    """walk the file once and return list of (header, aligned_seq) tuples.

    each entry begins with a '>' line; subsequent lines until the next '>'
    are concatenated (whitespace stripped) into the alignment string.
    """
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
                # strip any whitespace inside line (defensive)
                chunks.append(line.strip())
        # flush trailing record
        if header is not None:
            records.append((header, "".join(chunks)))
    return records


# ---------------------------------------------------------------------------
# residue-to-column mapping
# ---------------------------------------------------------------------------

def buildPositionToColumnMap(aligned_seq, target_positions):
    """walk aligned_seq, count non-gap chars, record column index for each
    requested 1-based residue position.

    returns dict {residue_position: column_index_0based}
    """
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


# ---------------------------------------------------------------------------
# tally amino acids at a column across all rows
# ---------------------------------------------------------------------------

def tallyColumn(records, col_idx):
    """count amino acids at col_idx across every record.

    returns (Counter excluding '-', n_gaps, n_total)
    n_total = number of sequences (some may be shorter than col_idx, in which
    case treat as 'X' = missing; we count them as gaps to be conservative).
    """
    counts = Counter()
    n_gaps = 0
    n_total = len(records)
    for _, seq in records:
        if col_idx >= len(seq):
            # sequence shorter than expected — treat as gap
            n_gaps += 1
            continue
        aa = seq[col_idx]
        if aa == "-":
            n_gaps += 1
        else:
            counts[aa] += 1
    return counts, n_gaps, n_total


# ---------------------------------------------------------------------------
# verdict thresholds
# ---------------------------------------------------------------------------

def conservationVerdict(frac_wt):
    if frac_wt >= 0.90:
        return "HIGHLY CONSERVED"
    if frac_wt >= 0.50:
        return "MODERATELY CONSERVED"
    if frac_wt >= 0.30:
        return "WEAKLY CONSERVED"
    return "VARIABLE"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print(f"[info] parsing alignment: {ALIGNMENT_PATH}")
    records = parseAlignedFasta(ALIGNMENT_PATH)
    n_seqs = len(records)
    print(f"[info] parsed {n_seqs} sequences")

    # alignment columns — use first record's length, then sanity-check uniformity
    aln_len = len(records[0][1])
    n_matching_aln_len = sum(1 for _, s in records if len(s) == aln_len)
    print(f"[info] alignment length (cols) = {aln_len}")
    print(f"[info] sequences matching alignment length = {n_matching_aln_len}/{n_seqs}")

    # locate EcDrt3b
    ref_idx = None
    for i, (hdr, _) in enumerate(records):
        if REF_ACCESSION in hdr:
            ref_idx = i
            break
    if ref_idx is None:
        raise SystemExit(f"[fatal] reference {REF_ACCESSION} not found in alignment")
    ref_header, ref_aligned = records[ref_idx]
    print(f"[info] EcDrt3b header: {ref_header}")

    # sanity check: gap-stripped reference
    ref_ungapped = ref_aligned.replace("-", "")
    print(f"[info] EcDrt3b gap-stripped length = {len(ref_ungapped)}")
    print(f"[info] EcDrt3b first 30 residues: {ref_ungapped[:30]}")
    if len(ref_ungapped) != REF_EXPECTED_LEN:
        print(f"[warn] expected length {REF_EXPECTED_LEN}, got {len(ref_ungapped)}")
    if not ref_ungapped.startswith(REF_EXPECTED_PREFIX):
        print(f"[warn] expected prefix {REF_EXPECTED_PREFIX!r}, got {ref_ungapped[:len(REF_EXPECTED_PREFIX)]!r}")

    # build position -> column mapping
    target_positions = [pos for pos, _, _ in TARGETS]
    pos_to_col = buildPositionToColumnMap(ref_aligned, target_positions)
    missing = [p for p in target_positions if p not in pos_to_col]
    if missing:
        raise SystemExit(f"[fatal] could not map residues: {missing}")

    # cross-check that the ref aa at each mapped column equals the expected wt aa
    print("\n[info] residue-to-column mapping (verifying reference aa):")
    for pos, expected_aa, role in TARGETS:
        col = pos_to_col[pos]
        observed = ref_aligned[col]
        flag = "OK" if observed == expected_aa else f"MISMATCH (got {observed})"
        print(f"  pos {pos:>4} -> col {col:>5}  expected {expected_aa}  observed {observed}  [{flag}]")

    # tally and prepare rows
    rows = []
    for pos, expected_aa, role in TARGETS:
        col = pos_to_col[pos]
        counts, n_gaps, n_total = tallyColumn(records, col)
        denom = n_total - n_gaps  # exclude gaps from frac_wt denominator
        wt_count = counts.get(expected_aa, 0)
        frac_wt = (wt_count / denom) if denom > 0 else 0.0
        gap_frac = (n_gaps / n_total) if n_total > 0 else 0.0
        top5 = counts.most_common(5)
        rows.append({
            "position": pos,
            "wt_aa": expected_aa,
            "role": role,
            "alignment_column": col,
            "n_total": n_total,
            "n_with_aa": wt_count,
            "n_gaps": n_gaps,
            "frac_wt": frac_wt,
            "gap_frac": gap_frac,
            "top5": top5,
            "verdict": conservationVerdict(frac_wt),
        })

    # stdout table
    print("\n" + "=" * 100)
    print("CONSERVATION TABLE (frac_wt excludes gaps from denominator)")
    print("=" * 100)
    header_fmt = "{:>4}  {:>5}  {:>6}  {:>7}  {:>9}  {:>6}  {:>7}  {}"
    row_fmt    = "{:>4}  {:>5}  {:>6}  {:>7}  {:>9}  {:>6}  {:>7.3f}  {}"
    print(header_fmt.format("Pos", "WT", "AlnCol", "N_total", "N_with_AA", "N_gaps", "Frac_WT", "Top5"))
    print("-" * 100)
    for r in rows:
        top5_str = ", ".join(f"{aa}:{n}" for aa, n in r["top5"])
        print(row_fmt.format(
            r["position"], r["wt_aa"], r["alignment_column"],
            r["n_total"], r["n_with_aa"], r["n_gaps"],
            r["frac_wt"], top5_str
        ))

    # CSV
    os.makedirs(os.path.dirname(CSV_OUT), exist_ok=True)
    csv_header = [
        "position", "wt_aa", "alignment_column",
        "n_total_seqs", "n_seqs_with_aa", "n_gaps",
        "frac_wt", "gap_frac",
        "top1_aa", "top1_count", "top1_frac",
        "top2_aa", "top2_count", "top2_frac",
        "top3_aa", "top3_count", "top3_frac",
        "top4_aa", "top4_count", "top4_frac",
        "top5_aa", "top5_count", "top5_frac",
        "verdict", "role",
    ]
    with open(CSV_OUT, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(csv_header)
        for r in rows:
            denom = r["n_total"] - r["n_gaps"]
            row = [
                r["position"], r["wt_aa"], r["alignment_column"],
                r["n_total"], r["n_with_aa"], r["n_gaps"],
                f"{r['frac_wt']:.6f}", f"{r['gap_frac']:.6f}",
            ]
            for i in range(5):
                if i < len(r["top5"]):
                    aa, n = r["top5"][i]
                    frac = (n / denom) if denom > 0 else 0.0
                    row.extend([aa, n, f"{frac:.6f}"])
                else:
                    row.extend(["", "", ""])
            row.extend([r["verdict"], r["role"]])
            writer.writerow(row)
    print(f"\n[info] CSV written to {CSV_OUT}")

    # verdicts
    print("\n" + "=" * 100)
    print("VERDICTS (frac_wt thresholds: >=.90 HIGHLY, .50-.90 MODERATELY, .30-.50 WEAKLY, <.30 VARIABLE)")
    print("=" * 100)
    for r in rows:
        print(
            f"  {r['wt_aa']}{r['position']:<4}  frac_wt={r['frac_wt']:.3f}  "
            f"gap={r['gap_frac']:.3f}  -> {r['verdict']}  ({r['role']})"
        )

    # summary
    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)
    print(f"  total sequences parsed       : {n_seqs}")
    print(f"  alignment columns            : {aln_len}")
    print(f"  seqs matching alignment len  : {n_matching_aln_len}")
    print(f"  reference                    : {REF_ACCESSION} (EcDrt3b)")
    print(f"  reference ungapped length    : {len(ref_ungapped)} (expected {REF_EXPECTED_LEN})")


if __name__ == "__main__":
    main()
