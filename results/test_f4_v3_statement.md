# Test F4 — Mode 2 simulation: drop-in v3 statement

## 1. L vs I_struct scaling for Mode 2

| L (codons) | I_struct (bits) | I per codon | H(Y) per codon |
|------------|------------------|-------------|------------------|
|    4 |       16.8632 |    4.2158 |        4.2191 |
|    8 |       33.7315 |    4.2164 |        4.2194 |
|   16 |       67.3605 |    4.2100 |        4.2135 |
|   32 |      134.8820 |    4.2151 |        4.2185 |
|   64 |      269.6645 |    4.2135 |        4.2171 |
|  128 |      539.5858 |    4.2155 |        4.2191 |
|  256 |     1078.8660 |    4.2143 |        4.2181 |
|  500 |     2107.0914 |    4.2142 |        4.2178 |

At eps = 1e-4 with the standard genetic code, I_struct per codon
saturates at the codon-count entropy H(Y) ≈ 4.2181 bits (closed form).
The 4.32-bits/codon log2(20) ceiling cited in the dispatch is the
*alphabet ceiling*; the *information* ceiling is the strictly-tighter
codon-count entropy, because the 64 codons are unevenly distributed
across the 21 amino-acid symbols.

## 2. Mode 1 vs Mode 2 capacity comparison at L_codons = 64

| L_codons | I_Mode1_DNA (3·L nt) | I_Mode2 peptide | ratio |
|----------|----------------------|-----------------|-------|
|   4 |           23.9813 |        16.8881 | 1.4200 |
|   8 |           47.9466 |        33.7194 | 1.4219 |
|  16 |           95.9015 |        67.4186 | 1.4225 |
|  32 |          191.8566 |       134.8146 | 1.4231 |
|  64 |          383.6882 |       269.7221 | 1.4225 |
| 128 |          767.4374 |       539.6943 | 1.4220 |
| 256 |         1534.7789 |      1079.2931 | 1.4220 |

Predicted ratio (closed form):  1.4219
Dispatch-cited ratio 6 / 4.32:  1.3883
Closed-form ratio 6 / 4.218:    1.4224

The empirical Mode-1-over-Mode-2 ratio matches the codon-count-entropy
prediction at L_codons ≥ 16, confirming P_F4_1: code degeneracy
strictly reduces information capacity from the substrate ceiling.

## 3. Apparatus signature for Mode 2 at L = 64

  I_struct^pop      = 269.6739 bits  (predicted ≈ 64 × 4.218 = 270 bits)
  I_struct^chan     = 0.5141 bits  (KL(Mode2_standard || ensemble mixture), L_chan=2, n=8000)
  mechanism fidelity = 0.999881        (predicted 1 - eps = 0.9999)

This confirms P_F4_3: the dual-observable apparatus extends cleanly
to the cross-descriptor templating case. Mode 2 is distinguishable
from Mode 1 DNA (different alphabets), from Mode 2 swapped (same
alphabet but different codon→aa mapping), and from AbiK uniform
(structured non-uniform marginal P(Y)).

## 4. Nested-Mode-1 copyability outcome

K = 200 agents, 500 generations, μ_meta = 1e-4 per gene bit.
Initial gene = uniform random 320-bit vector (random codon→aa table).
Selection on phenotype peptide vs fixed target.

Final-generation values (gen 500):
  mean fitness         = 0.4325
  max fitness          = 0.4333
  max achievable       = 0.8333
  gene diversity       = 0.0134
  code fidelity        = 0.5200

The population evolves toward a fitness-favorable code: phenotype
fitness rises and gene diversity drops. Code fidelity (fraction of
input-appearing codons whose population-modal aa matches the optimal
aa) rises from chance-level to a substantial fraction by generation
500. This is the operational test of P_F4_4: the meta-level code
*itself* is heritable via the Mode-1 inheritance of the gene.

## 5. Code expansion (20 → 24 aa) result

Standard code H(Y) per codon  = 4.2171 bits  (predicted 4.218)
Expanded 24-aa H(Y) per codon = computed at L=64 (see Step 5 csv)
  empirical ratio             = 1.0439
  closed-form ratio           = 1.0433
  naive log2(24)/log2(20)     = 1.0609

The expanded-alphabet ratio (1.044 empirical / 1.043 closed-form)
sits below the loose log₂(24)/log₂(20) ≈ 1.061 alphabet-ratio
ceiling. Reason: the dispatch's 1.061 prediction assumes the codons
redistribute uniformly across the 24 aa, but in the constructed
expansion only 3 codons are reassigned (one each from the top-3
over-represented aa: Leu, Ser, Arg → new labels), so the codon-count
distribution remains structured. The codon-count entropy of the
realized expanded mapping is 4.401 bits/codon, giving a tighter ratio
of 1.043. The increase is heritable only via Mode-1 inheritance of
the expanded charging machinery, confirming the qualitative claim of
P_F4_5 — code expansion increases I_struct, gated by Mode 1.

## 6. Drop-in paragraph for v3 Results section "Mode 2 (translation)"

> Mode 2 (code-encoded cross-descriptor templating, exemplified by
> ribosomal translation) is a distinct substrate in our classification:
> a 1D sequence template X over a 4-nucleotide alphabet is mapped
> through a separable, evolved lookup table T into a product Y over
> a 21-symbol amino-acid alphabet. We simulated Mode 2 explicitly at
> per-codon translation error ε = 10⁻⁴ and verified two quantitative
> claims. First, Mode 2's information capacity per template codon
> saturates at the codon-count entropy of T — for the standard genetic
> code this is H(Y) = 4.218 bits per codon, strictly less than the
> Mode-1 substrate ceiling of 6 bits per codon (3 nucleotides × 2
> bits). The signature degeneracy fingerprint (Leu/Ser/Arg over-
> represented, Met/Trp under-represented) is recovered at L = 64.
> Second, Mode 2's joint G2 apparatus signature (I_struct^pop,
> I_struct^chan, per-codon fidelity) cleanly separates it from
> Mode 1 DNA, from a permuted-code Mode 2, and from a uniform-output
> non-templating channel — confirming that the apparatus extends to
> cross-descriptor templating. Critically, in a two-level population-
> dynamics simulation where the lookup table T is itself a Mode-1-
> inherited gene (320 bits, μ = 10⁻⁴ per bit, K = 200, 500
> generations), selection on the *phenotype* peptide produces
> heritable convergence of the *meta-level code* toward a fitness-
> favorable mapping. Mode 2 inherits Mode 1's open-ended copyability
> via the gene encoding the operator, satisfying our inheritance
> theorem. An expanded 24-aa code (constructed by reassigning three
> codons from the most-degenerate standard aa to three new labels)
> yields I_struct per codon ~4.4% larger than the standard 21-symbol
> code (empirical 1.044, closed-form 1.043). This is below the loose
> alphabet-ratio ceiling log₂(24)/log₂(20) ≈ 1.061 because the
> realized codon-count distribution sets the tighter information
> ceiling. The heritability of any code expansion is gated by
> Mode-1 inheritance of the gene encoding the new charging machinery.

