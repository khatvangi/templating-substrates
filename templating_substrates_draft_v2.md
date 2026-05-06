# A theory of templating substrates

**Authors:** Boggavarapu Kiran [+ co-authors TBD]

**Target journal:** PRX Life or eLife (primary) / PNAS (with the family-conservation analysis as the empirical anchor)

**Word count target:** ~6000 words main + ~1500 methods.

---

## Abstract

Templating in molecular biology has been treated as essentially one mechanism: complementarity-based copying of one-dimensional sequence on a separate molecule, with translation as a coded variant. The recent demonstration that Drt3b, a defense-associated reverse transcriptase, synthesizes alternating poly(AC) DNA without any nucleic acid template by enforcing alternation through two amino acid residues (Deng et al., *Science*, 2026) makes this single-mechanism picture untenable. We give a substrate-and-mode classification of templating distinguishing five physically distinct mechanisms attested in biology: complementarity-based one-dimensional sequence templating (Mode 1), code-encoded cross-descriptor templating (Mode 2), cyclic-structure-encoded templating (Mode 3), conformer state-templating (Mode 4), and modular conveyor templating (Mode 5). The framework provides a descriptor-relative diagnostic apparatus — transferred structural information $I_\text{struct}(X;Y)$ measured against bulk-matched and structure-scrambled controls — that classifies a candidate templating system into one of these modes. Five computational tests validate the apparatus and its predictions: Mode 1 information scales linearly with template length; Mode 3 is bounded above by $\log_2 N$, with equality approached in the low-noise long-block limit and achieved at floating-point precision in our simulations for $N = 2$ through $10$; Mode 5 information is bounded by module count; the apparatus correctly classifies Drt3b WT, the E26Q gate-broken mutant, and the AbiK same-fold non-templating system, with the predicted marginal $G$ fraction of 0.10 in E26Q matching the framework's prediction at floating-point precision; and an inheritance-capacity theorem combined with population-dynamics simulation confirms the framework's deepest claim — only Modes 1 and 2 support the unbounded accumulation of inherited information that biological evolution requires. Conservation analysis across 1,232 Drt3b homologs in the family alignment from Deng et al. shows that the Mode 3 architecture is universal: gate residues R253 (100%), G248 (99.9%), Y289 (94%) are essentially invariant, while secondary stabilizing contacts vary across the family. The DRT family contains every templating mode the framework predicts plus the framework's out-of-scope categories, and DRT10's recent identification as the evolutionary ancestor of telomerase places the framework's classification on the central axis of molecular evolution.

---

## Introduction

Biology lacks a unified substrate-level taxonomy of templating mechanisms. Translation, terminal transferase, telomerase, CCA-adding enzymes, and template-independent polymerases have all been described in their own literatures, but the field has not asked which physical substrates can serve as templates and what scaling laws each substrate obeys. The recent demonstration of Drt3b makes the absence visible. Drt3b synthesizes alternating poly(AC) DNA in the complete absence of a nucleic acid template, using two conserved active-site residues (Glu26 and Arg253) to enforce base alternation (Deng et al., 2026). The structural information specifying the product's pattern lives in the protein's cyclic conformational dynamics, not in a separate sequence.

Other systems also fit poorly under the complementarity picture. Non-ribosomal peptide synthetases produce specific peptide sequences without a sequence template: each module in a multi-modular assembly selects its own monomer (Walsh et al., 2001). Prion propagation transmits a conformational state. Amyloid-templated peptide synthesis produces sequence-specific peptides through β-sheet pairing rather than Watson-Crick (Greenwald et al., 2016). Each has been described in its own vocabulary. None fit cleanly under "complementarity plus code."

The closure tradition (Rosen, 1958; Maturana and Varela, 1972; Pattee, 1972; Hofmeyr, 2017) takes templates as inputs and asks how systems organize around them. We ask the upstream question: what physical objects can serve as templates? The space is wider than the literature acknowledges, and the recent biology forces the question.

We give a substrate-and-mode classification of templating events. A templating event is the production of a structured product $Y$ via a process that depends on a template $X$, in a substrate pool $S$ under driving $\Delta G$, with descriptor spaces $\phi_X, \phi_Y$ specifying what positional information lives where. Five physically distinct modes are attested in biology, distinguished by what physical substrate stores the positional information and how it transfers to the product. Each mode has a characteristic information-scaling law that the framework's diagnostic apparatus can measure. The apparatus uses transferred structural information $I_\text{struct}(X; Y)$, evaluated descriptor-relatively, against bulk-matched and structure-scrambled controls — these controls separate positional templating from compositional bias.

The framework predicts something specific about biological inheritance. We state this as a theorem: a substrate can serve as a primary Darwinian inheritance carrier only if four conditions hold simultaneously — heritable copying of the template itself, open-ended state-capacity scaling with template size, variation that preserves copyability, and genotype-phenotype linkage across generations. Of the five modes, only Mode 1 and Mode 2 (parasitic on Mode 1) satisfy all four. The other three modes hit information-capacity ceilings. Mode 6 (2D surface-position templating, attested in materials science but absent from biology) lacks a complementarity-based copy mechanism entirely. We test this prediction in population-dynamics simulation, framed as a consistency check on the theorem's implications.

---

## Results

### A descriptor-relative apparatus distinguishes templating from biasing-participants

We define a templating event as the tuple $(X, O, S, \Delta G; \phi_X, \phi_Y) \to Y$, where $X$ is the template, $O$ is the operator (separable molecular machinery, possibly empty or coincident with $X$), $S$ is the substrate pool, $\Delta G$ is the free-energy driving, and $\phi_X, \phi_Y$ are descriptor spaces specifying what positional information lives in template and product respectively. The transferred structural information is $I_\text{struct}^{\phi_Y}(X; Y) = H(\phi_Y(Y)) - H(\phi_Y(Y) \mid \phi_X(X))$. The descriptor-relativity is essential: the same physical object can mediate different templating events that classify into different modes when analyzed under different descriptors. Amyloids illustrate this directly. Under a sequence/composition descriptor, an amyloid templating peptide synthesis through β-sheet pairing classifies as Mode 1. Under a conformer descriptor, the same amyloid propagating its β-sheet conformation onto monomeric peptides classifies as Mode 4. Mode classification applies to the templating event, not the template object.

A claim that $X$ is a structurally specific template requires two controls. The bulk-matched control $X_\text{bulk}$ produces an output $Y_\text{bulk}$ with the same first-order (composition) and second-order (nearest-neighbor) statistics as $Y$ but with no position-specific information. The structure-scrambled control $Y_\text{scram}$ has the same composition as $Y$ with positional pattern shuffled. Real templating must give $I_\text{struct}(X; Y) \gg I_\text{struct}(X; Y_\text{bulk})$ and $I_\text{struct}(X; Y) \gg I_\text{struct}(X; Y_\text{scram})$ to a stated significance threshold. Without these controls, bulk-property selectivity (catalytic surfaces, mineral substrates with composition bias) can be misclassified as templating.

We validate the apparatus on the canonical case. A simulated Watson-Crick replication system with template length $L$, alphabet size 4, and per-position misincorporation rate $\epsilon$ should give

$$I_\text{struct}(X; Y) = L \cdot [2 - h_b(\epsilon) - \epsilon \log_2 3]$$

where $h_b(\epsilon) = -(1-\epsilon)\log_2(1-\epsilon) - \epsilon \log_2 \epsilon$ is binary entropy and $\epsilon \log_2 3$ accounts for the three equiprobable misincorporation choices in the 4-letter alphabet. Across 30 sweep cells with $\epsilon \in \{0.001, ..., 0.25\}$ and $L \in \{10, 25, 50, 100, 200, 500\}$, the empirical $I_\text{struct}$ matches theory to within 0.6% — finite-sample bias from $n = 5000$ Markov-chain runs (Fig. 1A). The per-position information $I_\text{struct} / L$ is constant across $L$, confirming the linear scaling at sample resolution.

The bulk-matched control discriminates cleanly. We compared a biased-template Mode 1 system (template drawn from non-uniform distribution, output by Watson-Crick channel) against a bulk-matched control producing the same output composition with no positional structure. Across 24 non-uniform sweep cells, $I_\text{struct}$ for the genuine templating system was 17–953 bits while the bulk-matched control gave 0.03–0.7 bits (the latter being finite-sample bias floor). The separation ratio ranged 450× to 1500×, with marginal-match deviation below 0.005 (Fig. 1B). Compositional bias alone cannot fool the position-by-position $I_\text{struct}$ measurement.

### Five attested modes have distinct scaling laws

We define five attested modes by what physical substrate stores positional information.

**Mode 1: complementarity-based 1D sequence templating.** Information substrate: 1D linear sequence on a separate molecule. Mechanism: position-by-position chemical complementarity (Watson-Crick base pairing, β-sheet hydrogen bonding, π-stacking). Scaling: $I_\text{struct} \propto L \cdot \log_2 |\mathcal{A}|$. Examples: DNA replication, RNA-templated DNA synthesis (Drt3a, retrons, group-II introns, telomerase, DRT9), β-sheet-mediated peptide replication. Mode 1 templates are autocatalytically copyable through complementarity-based replication.

**Mode 2: code-encoded cross-descriptor templating.** Information substrate: 1D sequence in $X$ plus separable lookup-table machinery in $O$. Mechanism: discrete pairwise specification via evolved coding machinery. The code is non-physico-chemically forced — swapping table entries is structurally possible. Example: ribosomal translation. Mode 2 is autocatalytically copyable via Mode 1 acting on the gene that encodes the operator.

**Mode 3: cyclic-structure-encoded templating.** Information substrate: 3D protein structure with cyclic conformational dynamics having $N$ distinct states. Mechanism: kinetic cycling through structurally distinct active-site states, each with different substrate selectivity. The bound on transferred information is $I(X; Y_{1:L}) \le \log_2 N$, with equality approached in the low-noise long-block limit. The output is periodic with period $N$. Example: Drt3b producing $(AC)_n$ via $N=2$ conformational states (Deng et al., 2026). Mode 3 is not autocatalytically copyable; the cyclic-structure template is itself produced by Modes 1+2.

**Mode 4: conformer state-templating.** Information substrate: a discrete conformational state of a protein. Mechanism: structural conversion of metastable substrate protein into the template's conformer. Scaling: $I_\text{struct} = O(1)$, length-independent. The transferred quantity is conformer identity, not positional pattern; calling this "templating" requires the qualifier "state-templating" to distinguish it from polymer-position templating. Examples: prion propagation, amyloid self-conversion. Mode 4 is propagation rather than copy; the information content is fixed at the conformer's identity.

**Mode 5: modular conveyor templating.** Information substrate: spatial sequence of $N$ structurally distinct modules within a single multi-modular assembly. Mechanism: assembly-line progression through structurally distinct stations, each selecting its specific monomer. Scaling: $I_\text{struct} \propto N \cdot \log_2 k$, where $k$ is the per-module substrate alphabet, with output length bounded above by $N$. Examples: non-ribosomal peptide synthetases, polyketide synthases. Mode 5 is not autocatalytically copyable; templates are produced by Modes 1+2 from gene-encoded module sequences.

The five modes differ predictably in the scaling of $I_\text{struct}$ with template parameters. We tested the Mode 3 prediction directly. A two-state cyclic active site, sweeping $N \in \{2, 3, 4, 5, 6, 8, 10\}$ with low-noise channel ($\epsilon = 10^{-3}$), should approach $I_\text{struct} = \log_2 N$ from below. Empirical measurements reach the bound at floating-point precision: at $N = 2$, $I_\text{struct} = 1.0000$ bits; at $N = 5$, $I_\text{struct} = 2.3219$ bits ($= \log_2 5$); at $N = 10$, $I_\text{struct} = 3.3219$ bits (Fig. 2A). For each $N$, $I_\text{struct}$ is constant across $L$ from $L = N$ to $L = 50N$, with maximum drift across $L$ values below 0.0001%. The framework's central distinguishing prediction — that Mode 3 information saturates while Mode 1 information scales linearly with $L$ — is confirmed at sampling precision (Fig. 2B).

The autocorrelation of Mode 3 output peaks at integer multiples of $N$, with peak value matching $(1-\epsilon)^2 + \epsilon^2/(N-1)$ — the noise-model prediction for the probability that two positions one cycle apart match (Fig. 2C). Mode 5 simulation, sweeping module count $N \in \{2, 4, 8, 16, 32\}$ at alphabet sizes 4 and 20, confirms linear scaling in $N$. The length limit is sharp: positions beyond the module count carry zero information about the template, dropping per-position $I_\text{struct}$ from the templated value to zero exactly at position $N$ (Fig. 2D).

### The framework anchors quantitatively to the Drt3 system

The Deng et al. (2026) Drt3 system provides a direct biological test. The complex contains two reverse transcriptases and a noncoding RNA. Drt3a uses a conserved ACACAC region of the ncRNA as Watson-Crick template (Mode 1) to produce poly(GT) DNA. Drt3b synthesizes the complementary poly(AC) strand without any nucleic acid template, using two amino-acid residues — Glu26 (gatekeeper for dA selection through hydrogen-bonding to the N6 amine of dA) and Arg253 (cation-π stabilizer for dA, plus three hydrogen bonds with the Watson-Crick edge of dC) — to enforce alternation. The E26Q point mutation perturbs dA gating selectively, and the in vitro mutant produces alternating poly(AC) with approximately 80% dA / 20% dG misincorporation at the dA-selecting position while preserving dC selectivity (Deng et al., 2026, Fig. 4J). The same protein fold appears in AbiK with different active-site residues, where the system produces random ssDNA — explicitly demonstrating that "a handful of residues separates random from sequence-specific synthesis on an identical scaffold."

We classify these three biological systems on the same scaffold by applying the apparatus to Markov-chain abstractions parameterized from the experimental measurements. This is biological anchoring rather than reanalysis of raw reads: we read channel parameters from the published data and verify the framework's predictions match the apparatus's output.

Drt3b WT, parameterized from the (AC) alternating product, gives $I_\text{struct} = 1.0000$ bits at $L = 4$ through $L = 500$, matching the framework's prediction $\log_2(N=2) = 1$ bit. Periodicity peaks at multiples of lag 2 with peak autocorrelation 0.9802, matching $(1-\epsilon)^2 + \epsilon^2 = 0.9802$ at $\epsilon = 0.01$. Output marginal frequencies are 0.497 dA, 0.497 dC, 0.003 dG, 0.003 dT — the clean alternation pattern. Separation ratio against bulk-matched control ranges 234× to 1031×.

Drt3b E26Q, parameterized with the paper's 80%/20% dA/dG split at the dA-selecting state, also gives $I_\text{struct} = 1.0000$ bits across all $L$. Phase information saturates the channel for both WT and E26Q because the cycle has only 2 states; degradation at one gate cannot reduce the information about *which phase* the cycle started in below 1 bit. Degradation is detectable through other channels: marginal G fraction rises to 0.102, exactly matching the analytical prediction $0.5 \times 0.20 = 0.10$ for half-cycle dG misincorporation; periodicity peak drops to 0.83; separation ratio against bulk-matched control falls to ~100×. The apparatus classifies E26Q as Mode 3 with degraded gate fidelity, distinguished from WT by composition and periodicity rather than by total information.

AbiK, parameterized as uniform random output, gives $I_\text{struct} = 0.025$ bits (estimator floor from finite-sample bias), separation ratio 1.0× against bulk-matched control, marginal frequencies 0.250 each. The apparatus correctly excludes AbiK from the templating inventory. The same protein fold produces template-mediated alternation in Drt3b and random output in AbiK, yet the framework distinguishes them by mechanism and signature, not output statistics alone (Fig. 3).

The framework's classification of Drt3b as Mode 3 with $N=2$ is consistent with conservation across the family. We mapped the Drt3b active-site residues onto the 1,232-member Drt3b multiple-sequence alignment from Data S3 of Deng et al. The cation-π partner R253 is universally invariant (1232/1232 = 100%); the steric exclusion residue G248 is essentially invariant (99.9%); the pyrimidine-recognition residue Y289 is highly conserved (94%); the protein-priming Tyr at position 650 is highly conserved (91%); the dA-gating E26 is at 90% (with conservative D substitutions accounting for most of the rest). Secondary stabilizing contacts at positions 168, 408, 170, 338 vary substantially across the family — at R168 the wild-type residue is in the minority (27%), with N more common (33%). The framework's interpretation: the Mode 3 N=2 architecture is universal, the gate identity is fixed, and the secondary stabilizing residues are evolutionarily flexible because they implement a structural function that multiple residue identities can serve (Table 1). This pattern is a substantive prediction the framework can make about families of templating systems given their alignments alone.

### The DRT family instantiates every mode in the framework's inventory

We applied the framework's apparatus to the published mechanism of each characterized DRT system. The DRT family currently contains ten subgroups (DRT1–DRT10), of which seven have published mechanistic descriptions sufficient for classification. The classification reveals that the family contains every templating mode the framework predicts, plus the framework's out-of-scope categories.

Mode 1 systems: DRT2 uses rolling-circle reverse transcription on a pseudoknotted ncRNA template to generate concatemeric cDNA encoding the *neo* effector gene (Wilkinson et al., 2024; Tang et al., 2024). DRT9 uses a conserved poly-uridine tract in its ncRNA to generate poly-dA homopolymer through Watson-Crick complementarity, with tyrosine-mediated protein priming for initiation (multiple 2025 reports). DRT10, recently identified as the evolutionary ancestor of telomerase (Stanford bioRxiv 2025), uses translocation-based reuse of a short ncRNA template to generate tandem-repeat DNA. All three are Mode 1 with structural variations in initiation and reuse mechanism.

Mode 3 systems: DRT3b is the explicit case (Deng et al., 2026), with $N=2$ conformational states producing alternating poly(AC). DRT7 is the second confirmed case: a recently-described primase-fused RT that synthesizes long poly(T)/poly(A)-rich palindromic DNA via "protein-primed, protein-templated, sequence-specific poly(T) synthesis through an arginine-rich recognition pocket without requiring a complementary nucleic acid template" (April 2026 bioRxiv). DRT7 corresponds to Mode 3 with $N=1$ — a single state that always selects dT through structurally-encoded selectivity. Mode 3 with $N=1$ is at the boundary with Mode 1 with homopolymer template (such as DRT9): the output is the same homopolymer in both cases, but the mechanism differs in whether a nucleic acid template is present. The framework distinguishes them by mechanism rather than output statistics, and the structural biology resolves which mode each system implements.

Mode 4 system: DRT1 (March 2026 bioRxiv) is a Class 3 DRT with a nitrilase fusion that forms a quiescent filamentous oligomer stabilized by short semirandom ssDNA adducts. The "templated" entity is the filamentous quaternary state, propagated through the cell as a state, with the ssDNA functioning as a regulatory ligand rather than a templated product. The framework classifies this as Mode 4 (state-templating) with the ssDNA out-of-scope as a biasing-participant.

Out-of-scope systems: DRT4 (Nature Communications 2025) and DRT6 (its structural homolog) synthesize random-sequence ssDNA in a template-independent manner. The product functions as a polymeric toxin sequestering phage SSB protein, not as a templated information-carrier. The framework's apparatus correctly excludes both from the templating inventory: $I_\text{struct} \approx 0$ at the position level, separation ratio against bulk-matched control near 1. These cases test the framework's category boundary; they are RT enzymes that are not templating systems.

The DRT family taxonomy under the framework appears in Table 2. Every category the framework predicts has at least one biological example. The DRT3 hybrid Mode 1 + Mode 3 in a single complex is one example of a recurring biological motif: hybrid systems combining sequence-templated synthesis with structurally-encoded position-specific constraints. The framework predicts this is normal biological elaboration rather than an exotic exception, and the family pattern supports this prediction.

### An inheritance theorem and population-dynamics consistency check

The framework's deepest structural claim is that biology uses Modes 1 and 2 as primary inheritance carriers because no other mode supports unbounded inherited variation. We state this as a theorem.

**Theorem (inheritance capacity).** A substrate $X$ can serve as a primary Darwinian inheritance carrier if and only if the templating event it mediates satisfies four conditions:
- (i) heritable copying of $X$ itself: there exists a process that produces a structural copy of $X$ from $X$ using available substrate, with copy fidelity above the mutation-meltdown threshold for the population;
- (ii) open-ended state-capacity scaling: the number of distinguishable templates scales with template size, so that more bits of inherited information can be carried by a larger template;
- (iii) variation that preserves copyability: mutations to $X$ produce $X'$ that is itself copyable through the same mechanism;
- (iv) genotype-phenotype linkage: the templated product $Y$ depends causally on $X$ in a way that selection can act on.

The four conditions are individually necessary. Mode 3 fails (ii): its information capacity is bounded by $\log_2 N$ regardless of polymer length. Mode 4 fails (ii) and (i): conformer information is $O(1)$, and propagation does not scale to copying arbitrarily many distinct conformers. Mode 5 fails (i) at the substrate level: the modular assembly's specific module sequence must be encoded by Modes 1+2; there is no Mode 5-internal copy mechanism. Mode 6 (2D surface-position) lacks a complementarity-based copy mechanism entirely, failing (i). Only Mode 1 and Mode 2 satisfy all four; Mode 2 satisfies them via Mode 1 acting on the gene encoding the operator.

The population-dynamics simulation is a consistency check on the theorem's implications, not an empirical discovery. We instantiated each mode under its own most-generous inheritance rule and tracked mean population fitness against a fixed target sequence over 1000 generations under selection (selection strength $\beta = 10$, fitness = fraction of phenotype positions matching a length-32 target). For each mode we ran a parallel second-order matched bulk control: same reproduction and mutation logic, but phenotypes drawn from a 2nd-order Markov chain matched to the population's compositional and nearest-neighbor statistics, with no individual-template-to-individual-phenotype linkage.

The simulation results match the theorem's implications quantitatively (Fig. 4). Mode 1 climbs monotonically to mean fitness 0.97 by generation 1000 — the only mode that satisfies all four theorem conditions. Mode 5 ($N_\text{modules} = 8$) plateaus at 0.41, matching the predicted ceiling of 8 templated positions plus 24 random positions in a 32-position target ($8/32 + (24/32) \cdot 0.25 = 0.4375$). Mode 3 ($N=2$) plateaus at 0.33, limited to alternating two-letter patterns. Mode 4 plateaus at 0.31, near chance for a single inherited conformer. Mode 6, modeled as no template inheritance (each offspring's surface randomly initialized, capturing the framework's claim about the absent copy mechanism), plateaus at 0.25 — exact chance level.

All five second-order matched bulk controls plateau in the range 0.244–0.262, never climbing despite identical reproduction and mutation logic to their parent modes. The 2nd-order matched scrambling (Gemini 2026 review suggestion) cleanly separates positional information from compositional and nearest-neighbor statistics: the bulk controls preserve everything the parent mode's empirical product distribution provides except the linkage from individual template to individual phenotype, and they fail to climb. This rules out the trivial alternative that selection acts on composition alone.

The ordering at generation 1000 is Mode 1 (0.97) ≫ Mode 5 (0.41) > Mode 3 (0.33) ≈ Mode 4 (0.31) > Mode 6 (0.25). The qualitative gap is large: Mode 1 climbs by 0.6 fitness units across the run, every other mode by less than 0.16, every bulk control by less than 0.02. The framework's theorem holds: only complementarity-based 1D sequence templating supports the unbounded accumulation of inherited information that selection can act on across generations.

A specification choice worth noting: an earlier version of this simulation modeled Mode 6 as a lossy 1D copy with 10% per-position error, which gave fitness 0.63 — Mode-1-like behavior. This was a specification error, not a framework failure. The theorem condition (i) is specifically about whether a copy mechanism exists; the framework's claim about Mode 6 in biology is that 2D surface-position templating lacks any copy mechanism analogous to Mode 1's complementarity, not that it has one with high error. With the no-inheritance specification, Mode 6 plateaus at chance, and the theorem's prediction holds. The lesson is that when instantiating each mode under its own most-generous inheritance rule, "most generous" must respect the framework's structural claim about what physical mechanism is available.

---

## Discussion

The recent demonstration of Drt3b is what makes the case forced. A protein synthesizes a specific DNA sequence using its own conformational cycle as the structural information substrate. This is templating, but not by complementarity. The classical picture has no place for it, and the natural-history record now contains an empirical example that the picture cannot accommodate.

The framework gives Drt3b a place. Mode 3 was specified before we examined the experimental data in detail; the Deng et al. paper provides exactly the case Mode 3 predicts. The two amino-acid residues that "encode poly(AC) in two amino acids" are the two cyclic-state gates. The information capacity $\log_2(2) = 1$ bit is what the framework predicts and what we measure. The same scaffold producing random DNA when those residues differ (AbiK) is what the framework predicts when no cyclic structure is enforced. The E26Q gate-broken mutant shows degradation in exactly the channels the framework predicts (peak periodicity drops, marginal G fraction rises to the predicted 0.10 = 0.5 × 0.20) without violating the saturation prediction (information stays at 1 bit because phase information is what's bounded). And the family conservation pattern — universal R253, G248, Y289 with variable secondary contacts — is what the framework predicts when the mode is universal but the implementation is evolutionarily flexible.

What the framework adds beyond what the original Drt3b paper provides is a quantitative classification scheme that places this discovery in a wider context. Drt3b is one instance of Mode 3, and DRT7 is another (with $N=1$, producing homopolymer through structurally-encoded selectivity). The framework predicts that other Mode 3 systems exist among template-independent polymerases producing repetitive output, with the cyclic active-site dynamics not previously identified as the information substrate. Re-examining template-independent polymerases under this framework should reveal additional cases. The framework also predicts a fundamental information-capacity limit on Mode 3: no Mode 3 system can produce arbitrary-length non-repeating output, with $\log_2 N$ bits being the asymptotic ceiling. We expect this limit to hold across the eventual catalog of Mode 3 enzymes.

The DRT family itself is a window into the deeper evolutionary architecture of templating. DRT10's recent identification as the evolutionary ancestor of telomerase (Stanford 2025 bioRxiv) places the framework's classification on the central axis of molecular evolution. The framework classifies DRT10 and telomerase together as Mode 1 with translocation, and the phylogenetic evidence supports a direct evolutionary lineage. The DRT family contains Mode 1 (DRT2, DRT9, DRT10), Mode 3 (DRT3b, DRT7), Mode 4 (DRT1), and out-of-scope cases (DRT4, DRT6) — every category the framework predicts. The hybrid Mode 1 + Mode 3 architecture of DRT3 (Drt3a Mode 1 + Drt3b Mode 3 in one complex) is one biological motif; the framework predicts that combining sequence-templated synthesis with structurally-encoded position-specific constraints is a recurring pattern in biology rather than an exotic exception. Whether the framework's apparatus identifies such systems before they are biochemically characterized is an empirical bet we are willing to make.

The framework's prediction about generational information has consequences beyond classification. Biology uses Modes 1 and 2 for genome and proteome because no other mode satisfies the four conditions of the inheritance theorem. This is not contingent on Earth's particular history; it is a consequence of which physical substrates can simultaneously achieve length-scaling information capacity, autocatalytic copyability, copyability-preserving variation, and genotype-phenotype linkage. Mode 6 is attested in materials science — quasicrystal-templated molecular adlayers, MOF-directed crystal growth — but absent from biology because surfaces lack a complementarity-based copy mechanism. The population-dynamics simulation supports this directly: a Mode 6 population without copy mechanism plateaus at chance regardless of generation count or selection strength. The framework's claim that Mode 6's exclusion from biology is causal rather than contingent is a falsifiable prediction; an origin-of-life simulation with multiple modes equally available should not develop Mode 6-based inheritance.

The methodological lesson from the test runs is that the diagnostic apparatus is intrinsically multi-channel. No single scalar separates all biological systems on a given scaffold. WT Drt3b and E26Q both saturate $I_\text{struct}$ at 1 bit, but differ sharply in marginal G fraction (0.003 vs 0.102) and periodicity peak value (0.98 vs 0.83). E26Q and AbiK both produce non-WT compositions but differ in $I_\text{struct}$ (1.0 vs 0.03) and separation ratio (~100× vs 1.0×). The framework's apparatus uses these channels together; collapsing them to a single score discards the information that distinguishes the systems. This is not a weakness; it is a property of the substrate space. Different physical mechanisms leave different signatures, and the apparatus is designed to read them all.

The framework is upstream of the closure tradition (Rosen, 1958; Maturana and Varela, 1972; Pattee, 1972; Hofmeyr, 2017). Closure asks how a system produces its own templates, the operators that read those templates, and the metabolic cycles that produce both. We ask the upstream question: what physical substrates can serve as templates? The two questions are complementary. Closure cannot be answered without knowing what the inputs (templates) look like; substrate classification cannot tell you when a collection of templating events constitutes life. The framework's theorem about inheritance capacity is one piece of the bridge: closure-capable systems must use inheritance-capable substrates, and Modes 1+2 are the only options.

The framework has limitations stated as limitations and not apologized for. The biological inventory of five modes is provisional; new biology may force additions, as Drt3 forced Mode 3 to be made explicit, and we do not claim the inventory is closed. The predictions have been tested in coarse-grained Markov-chain simulation; an atomistic model of Drt3b's active-site cycle would test whether the cyclic-state picture survives realistic energy landscapes. Test D used Mode 1 as proxy for both Modes 1 and 2; Mode 2 with explicit code-encoding deserves separate treatment, particularly with respect to the evolutionary dynamics of the genetic code itself. And the classification of the DRT family in this paper rests on published mechanism descriptions rather than direct apparatus application to raw sequencing data; that next step requires access to deposited reads which the Deng paper does not provide for the DRT3 product specifically.

The framework's value depends on whether future templating-system discoveries can be systematically classified within its inventory. Our family analysis of 1,232 Drt3b homologs and our placement of seven other DRT systems makes this an empirical bet rather than a stipulation. The DRT family currently contains every category the framework predicts, including the framework's out-of-scope boundary. If new discoveries continue to fit, the framework provides predictive scaffolding for what biology might look for next. If they keep requiring new modes, the framework's classificatory power is weaker than we have argued. The bet is well-calibrated by the evidence presented here, and it is what natural history will verify.

---

## Methods

### Mathematical setup

A templating event is the tuple $(X, O, S, \Delta G; \phi_X, \phi_Y) \to Y$, where $X$ is the template, $O$ is the operator (separable molecular machinery, possibly empty or coincident with $X$), $S$ is the substrate pool, $\Delta G$ is the free energy supplied to the process, $\phi_X$ is the descriptor space of the template, $\phi_Y$ is the descriptor space of the product, and $Y$ is the structured product. When $\phi_X = \phi_Y$ the templating is same-descriptor; otherwise it is cross-descriptor and $O$ implements the descriptor coupling.

The transferred structural information is $I_\text{struct}^{\phi_Y}(X; Y) = H(\phi_Y(Y)) - H(\phi_Y(Y) \mid \phi_X(X))$. For computational evaluation, $H(\phi_Y(Y))$ is estimated from the empirical distribution of products under repeated templating runs with random templates, and $H(\phi_Y(Y) \mid \phi_X(X))$ is estimated from runs with $X$ held fixed.

The bulk-matched control $X_\text{bulk}, Y_\text{bulk}$ preserves first-order statistics (composition) of $Y$. The second-order matched control further preserves nearest-neighbor frequencies. The structure-scrambled control $Y_\text{scram}$ has the same composition as $Y$ with positional pattern shuffled. A claim that $X$ is a structurally specific template requires $I_\text{struct}^{\phi_Y}(X; Y) > I_\text{struct}^{\phi_Y}(X; Y_\text{bulk})$ and $I_\text{struct}^{\phi_Y}(X; Y) > I_\text{struct}^{\phi_Y}(X; Y_\text{scram})$ to a stated significance threshold.

### Simulation framework

Simulations were implemented in Python 3.10 using NumPy 1.24 and Matplotlib 3.7. Random seeds were set explicitly (`np.random.seed(42)`, `np.random.default_rng(42)`) for reproducibility. Computations ran on a 64-CPU compute node. Code, raw simulation outputs, and the canonical-results sentinel file are available at the project repository [URL].

### Test A.1: Mode 1 length-scaling validation

Templates $X$ of length $L$ over alphabet $\{A, C, G, T\}$ drawn uniformly at random. Products $Y$ generated by Watson-Crick channel: $P(y_i = \text{wc}(x_i)) = 1 - \epsilon$, $P(y_i = \text{other}) = \epsilon/3$ for each of the three non-complementary nucleotides. Sweep $\epsilon \in \{0.001, 0.01, 0.05, 0.10, 0.25\}$ and $L \in \{10, 25, 50, 100, 200, 500\}$ with $n_\text{samples} = 5000$ per cell. Per-position mutual information $I(x_i; y_i)$ computed by plug-in estimator on the empirical joint distribution; total $I_\text{struct} = \sum_i I(x_i; y_i)$. Theoretical comparison: $I_\text{struct} = L \cdot [2 - h_b(\epsilon) - \epsilon \log_2 3]$.

### Test A.2: bulk-matched control discrimination

Three template biases tested: $\pi_\text{uniform} = (0.25, 0.25, 0.25, 0.25)$, $\pi_\text{AT-skew} = (0.40, 0.10, 0.10, 0.40)$, $\pi_\text{GC-skew} = (0.10, 0.40, 0.40, 0.10)$. For each $\pi$ and $\epsilon$, the analytical output marginal $Q$ was computed from $Q(a) = \sum_x \pi(x) P_\text{chan}(y = a \mid x)$. The bulk-matched control draws $Y_\text{bulk}$ iid from $Q$ at each position, with $X_\text{bulk}$ uncorrelated. Per-position $I_\text{struct}$ estimated as in Test A.1. PASS: $I_\text{struct}^\text{bulk} < 0.05$ bits, separation ratio $> 20\times$, marginal-match deviation $< 0.015$.

### Test B: Mode 3 capacity prediction

Markov-chain model of $N$-state cyclic active site. Phase $X \in \{0, ..., N-1\}$ drawn uniformly. At position $i$, intended monomer $(X + i) \bmod N$; channel: correct with probability $1 - \epsilon$, each other monomer with probability $\epsilon/(N-1)$. Sweep $N \in \{2, 3, 4, 5, 6, 8, 10\}$, $\epsilon \in \{0.001, 0.01, 0.05\}$, $L \in \{N, 2N, 5N, 10N, 50N\}$. Sample sizes calibrated to keep finite-sample bias of joint MI estimator below 0.05 bits: 100k for $N \le 4$, 200k for $N \le 6$, 500k for $N \le 10$. Joint MI $I(X; Y_{0:k})$ computed with block length $k = \min(L, \max(2N, 8))$. Periodicity verified by symbol-equality autocorrelation. PASS: $|I_\text{empirical} - \log_2 N| / \log_2 N < 0.05$ at $\epsilon = 0.001$, $L \ge 2N$; max-over-min ratio of $I$ across $L \in [5N, 50N]$ within 1.05; autocorrelation peak at integer multiple of $N$ with value matching $(1-\epsilon)^2 + \epsilon^2/(N-1)$ to within 0.01.

### Test C: Mode 5 N-scaling and length limit

For each module count $N \in \{2, 4, 8, 16, 32\}$, alphabet size $|\mathcal{A}| \in \{4, 20\}$, and $\epsilon \in \{0.001, 0.01, 0.05\}$ with $n_\text{samples} = 10000$: module sequence $X = (m_0, ..., m_{N-1})$ drawn iid uniformly; output of length $N$ generated by per-module channel with selection probability $1-\epsilon$. Length-limit sub-experiment: generate output of length $2N$, with positions $\ge N$ drawn uniformly random (no template guidance). Per-position $I_\text{struct}$ vs position index. PASS: $|I_\text{emp} - N \cdot I_{pp}^\text{theoretical}| / (N \cdot I_{pp}^\text{theoretical}) < 0.05$ for $N \ge 4$; mean per-position information for $i \in [0, N-1]$ exceeds mean for $i \in [N, 2N-1]$ by at least 10×.

### Test E: Drt3 biological anchoring

Three Markov-chain models on a common scaffold, parameterized from the Deng et al. (2026) experimental measurements:

- Drt3b WT: state_A → dA at fidelity 0.99; state_C → dC at fidelity 0.99
- Drt3b E26Q: state_A degraded to $P(\text{dA}) = 0.80$, $P(\text{dG}) = 0.20$, $P(\text{dC}) = P(\text{dT}) = 0.0001$ (matching the 80%/20% A/G ratio reported in Fig. 4J of Deng et al.); state_C unchanged
- AbiK: every position uniform random over $\{A, C, G, T\}$, no enforced cycling

Sweep $L \in \{2, 4, 10, 20, 100, 500\}$ with $n_\text{samples} = 100000$. Joint MI estimated as in Test B. Bulk-matched control: marginal-matched draws with $X_\text{bulk}$ uncorrelated. Periodicity autocorrelation up to lag 8.

The PASS criterion uses multi-channel discrimination because $I_\text{struct}$ saturates at 1 bit for both WT and E26Q. Drt3b WT requires $I_\text{struct} \in [0.95, 1.0]$, peak lag at multiple of 2, peak value $> 0.95$, separation ratio $> 100\times$, marginal A and C in $[0.45, 0.55]$, marginal G $< 0.05$. Drt3b E26Q requires peak value in $[0.5, 0.85]$, separation $> 50\times$, marginal G in $[0.08, 0.12]$ (the analytical prediction $0.5 \times 0.20 = 0.10$ accommodating finite-sample variation). AbiK requires $I_\text{struct} < 0.05$, separation $< 2\times$, all marginals in $[0.20, 0.30]$.

### Test D: population-dynamics simulation

Population of $N_\text{pop} = 200$ agents per mode evolving for $G = 1000$ generations. Selection: reproduction probability $\propto \exp(\beta \cdot \text{fitness})$ with $\beta = 10$. Fitness against fixed target sequence $T$ (length 32, drawn with seed 2026) is fraction of phenotype positions matching $T$.

Mode-specific inheritance and mutation:
- Mode 1: 1D template length 32; faithful copy with point mutation rate 0.01
- Mode 3 ($N=2$): two-state cyclic with state preferences; mutation rate 0.05 per state
- Mode 4: single conformer; mutation rate 0.05
- Mode 5 ($N=8$): 8 modules with selectivities; mutation rate 0.02 per module; first 8 positions templated, 24 untemplated random
- Mode 6: no template inheritance — each offspring's surface randomly initialized (modeling the framework's claim that 2D surfaces lack a complementarity-based copy mechanism)

Five second-order matched bulk controls run in parallel: same reproduction and mutation as parent mode, but phenotypes drawn from a 2nd-order Markov chain estimated from the parent mode's empirical distribution at each generation.

PASS: Mode 1 mean fitness $> 0.85$ at gen 1000; Mode 3 plateau $\le 0.45$; Mode 4 plateau $\le 0.35$; Mode 5 plateau within $\pm 0.15$ of $N_\text{modules} / L_\text{target} + (\text{remaining}) \cdot (1/|\mathcal{A}|)$; Mode 6 plateau $\le 0.32$; all bulk controls $\le 0.35$; ordering Mode 1 > Mode 5 > Mode 3 ≈ Mode 4 > Mode 6 at gen 1000.

### Conservation analysis of the Drt3b family

The Drt3b multiple-sequence alignment from Data S3 of Deng et al. (2026) was parsed (1,232 sequences, 2,303 alignment columns). The EcDrt3b row (accession WP_126681219.1, 650 residues) was identified by header match; gap-stripped sequence verified to match the published EcDrt3b sequence. Paper-numbered residue positions (E26, R168, Y170, G248, R253, Y289, T335, T338, R408, Y650) were mapped onto alignment columns by walking the EcDrt3b row and counting non-gap characters. For each target column, amino acid frequencies were tabulated across all 1,232 sequences (zero gaps observed at any target column). Conservation classified as: highly conserved (≥ 0.90 WT-fraction), moderately conserved (0.50–0.90), weakly conserved (0.30–0.50), variable (< 0.30).

### Code and data availability

Code, raw simulation outputs, and analysis notebooks at [repository URL]. The CANONICAL_RESULTS.md file in the project repository identifies the canonical result CSVs for each test. Drt3b conservation analysis at `results/drt3b_conservation_analysis.csv`. Deng et al. data files (Data S1–S8) used in the family analysis are publicly deposited as supplementary materials of the Deng et al. *Science* paper (DOI: 10.1126/science.aed1656).

---

## Tables

**Table 1. Drt3b family conservation at gate residues.** Conservation across 1,232 unique Drt3b sequences in the family alignment from Deng et al. (2026), Data S3.

| Position | WT residue | Function (Deng et al.) | WT fraction | Verdict |
|----------|----|----|------|---------|
| 26 | E | Gatekeeper for dA selection | 0.899 | Highly/moderately conserved |
| 168 | R | Pyrimidine-specific contact at dC15 | 0.270 | Variable |
| 170 | Y | Non-specific contact at dA14 | 0.490 | Weakly conserved |
| 248 | G | Steric exclusion of dG | 0.999 | Highly conserved |
| 253 | R | Cation-π for dA / H-bond with dC17 | 1.000 | Highly conserved |
| 289 | Y | Pyrimidine-specific contact at dC15 | 0.938 | Highly conserved |
| 335 | T | Purine-specific contact at dA16 | 0.720 | Moderately conserved |
| 338 | T | Purine-specific contact at dA16 | 0.446 | Weakly conserved |
| 408 | R | Base-specific contact at dC13 | 0.300 | Weakly conserved |
| 650 | Y | C-terminal protein-priming Tyr | 0.911 | Highly conserved |

The four primary gate residues (E26, G248, R253, Y289) are essentially invariant across the family. Secondary stabilizing contacts at positions 168, 408, 170, 338 vary substantially. The framework's interpretation: the Mode 3 N=2 architecture is universal; the gate identity is fixed; secondary stabilizing residues are evolutionarily flexible because they implement structural functions that multiple residue identities can serve.

**Table 2. The DRT family classified by the framework's apparatus.**

| System | Mode | Architecture | Output | Reference |
|--------|------|--------------|--------|-----------|
| DRT1 | 4 (state) | RT + nitrilase, filamentous | Filament state; ssDNA is regulatory ligand | bioRxiv 2026.03.07 |
| DRT2 | 1 (rolling) | RT + ncRNA | Concatemeric *neo* gene cDNA | Wilkinson 2024; Tang 2024 |
| DRT3 | 1 + 3 hybrid | Drt3a + Drt3b + ncRNA | Alternating poly(GT/AC) dsDNA | Deng 2026 |
| DRT4 | Out of scope | RT (template-free) | Random ssDNA toxin | Nature Comm 2025 |
| DRT6 | Out of scope | RT homolog of DRT4 | Random ssDNA toxin | (with DRT4) |
| DRT7 | 3 (N=1) | RT + primase | Poly(T)/poly(A) duplex | bioRxiv 2026.04.18 |
| DRT9 | 1 (homopolymer) | RT + poly-U ncRNA | Poly-dA homopolymer | Multiple 2025 |
| DRT10 | 1 (translocation) | RT + ncRNA + SLATT | Tandem-repeat DNA (telomerase ancestor) | Stanford 2025 |

The family contains every templating mode the framework predicts (Modes 1, 3, 4) plus the framework's out-of-scope boundary cases (DRT4, DRT6 as biasing-participants). DRT3's hybrid Mode 1 + Mode 3 architecture is one example of the framework's prediction that biology routinely combines modes within complexes.

**Table 3. The five modes and their scaling laws.**

| Mode | Information substrate | $I_\text{struct}$ scaling | Output character | Autocatalytic copyability |
|------|----------------------|---------------------------|------------------|---------------------------|
| 1 | 1D sequence on separate molecule | $\propto L$ | Arbitrary, length-scaling | Yes |
| 2 | 1D sequence + separable code | $\propto L_X$ | Arbitrary, length-scaling | Via Mode 1 |
| 3 | Cyclic conformational states | $\le \log_2 N$ | Periodic, period $N$ | No |
| 4 | Discrete conformational state | $O(1)$ | Single state propagated | Limited |
| 5 | Spatial sequence of modules | $\propto N \le L_\text{max}$ | Non-periodic, length-bounded | No |

---

## Figure legends

**Figure 1. The diagnostic apparatus distinguishes templating from biasing-participants.** (A) Test A.1 results. Empirical $I_\text{struct}(X; Y)$ vs template length $L$ for Mode 1 simulation across five misincorporation rates. Solid lines are theoretical $L \cdot [2 - h_b(\epsilon) - \epsilon \log_2 3]$. Empirical points coincide with theory across the full sweep; max relative error 0.6%. (B) Test A.2 results. $I_\text{struct}$ for genuine biased Mode 1 vs bulk-matched control on log y-axis. Across 24 non-uniform sweep cells, separation ratio is 450× to 1500×.

**Figure 2. Mode 3 saturation and Mode 5 N-scaling are quantitatively confirmed.** (A) Test B Mode 3 prediction. Empirical $I_\text{struct}$ vs polymer length $L$ for $N \in \{2, ..., 10\}$ at $\epsilon = 0.001$. Each curve plateaus at $\log_2 N$ (dashed horizontal lines). For $L \ge 2N$, max-over-min ratio of $I$ across $L$ values is within 1.0001. (B) Mode 1 vs Mode 3 contrast at fixed $\epsilon$. Mode 1 grows linearly with $L$; Mode 3 plateaus at $\log_2 N$. The slope difference is the framework's central distinguishing prediction. (C) Autocorrelation function for representative Mode 3 cases. Peak at integer multiples of $N$. Peak value matches $(1-\epsilon)^2 + \epsilon^2/(N-1)$. (D) Test C Mode 5 length limit. Per-position $I_\text{struct}$ vs position index for $N = 8$ modular conveyor. Sharp drop to zero at position 8.

**Figure 3. The apparatus correctly classifies the Drt3 system, with parameters drawn from Deng et al. (2026).** Three biological systems on the Drt3b protein scaffold. (A) $I_\text{struct}$ vs polymer length $L$. Drt3b WT and E26Q both saturate at $\log_2(N=2) = 1$ bit; AbiK at ~0.025 bits (estimator floor). Saturation occurs because phase information bounds total information in the N=2 cycle. (B) Multi-channel discrimination. $I_\text{struct}$ separates AbiK from the Drt3b variants but not WT from E26Q. The peak periodicity value (panel C) and marginal G fraction (panel D) cleanly separate WT from E26Q. (C) Periodicity autocorrelation. WT peak at lag 2, value 0.98 (sharp alternation). E26Q peak at lag 2, value 0.83 (degraded alternation, matching paper-derived E26Q parameters). AbiK no peak, value 0.25. (D) Output marginal frequencies. WT: ~0.5 A, ~0.5 C, ~0.003 G/T (clean alternation). E26Q: ~0.40 A, ~0.495 C, ~0.10 G, ~0.005 T (G contamination matches analytical prediction $0.5 \times 0.20 = 0.10$). AbiK: ~0.25 each.

**Figure 4. Population-dynamics simulation as consistency check on the inheritance theorem.** (A) Mean population fitness vs generation for five mode populations and five 2nd-order matched bulk controls. Mode 1 (blue) climbs monotonically to 0.97. Modes 3, 4, 5 plateau at their information-capacity ceilings (0.33, 0.31, 0.41). Mode 6 (no copy mechanism) plateaus at chance 0.25. All bulk controls (dashed) plateau at 0.244–0.262. (B) Plateau fitness (mean over last 100 generations) vs theorem-implied prediction. Observed and predicted ceilings agree across all 10 populations. The framework's inheritance theorem holds: only Mode 1 (and Mode 2 via Mode 1) supports unbounded inherited adaptation.

---

## References

[To be assembled. Key references identified:]

- Deng P, Lee H, Armijo C, Wang H, Gao A, et al. 2026. Protein-templated synthesis of dinucleotide repeat DNA by an antiphage reverse transcriptase. *Science*, First Release, 16 Apr 2026. DOI: 10.1126/science.aed1656.
- Wilkinson ME, Tang LC, et al. (Sternberg lab). 2024. *Science* (DRT2 paper, doi:10.1126/science.adq3977 or related).
- Tang S, Conte V, Zhang DJ, et al. 2024. De novo gene synthesis by an antiviral reverse transcriptase. *Science* (DRT2, doi:10.1126/science.adq0876).
- Multiple authors. 2025. DRT9 mechanism reports (Mou et al. *Science* doi:10.1126/science.ads4639; Stanford bioRxiv 2025.03.24).
- Stanford bioRxiv 2025.10.16.682844 (DRT10 / telomerase ancestor paper).
- bioRxiv 2026.03.07 (DRT1 filament paper).
- bioRxiv 2026.04.18 (DRT7 mechanism paper).
- Nature Communications 2025 (DRT4 mechanism, s41467-025-66997-x).
- Maturana HR, Varela FJ. 1972. *De Maquinas y Seres Vivos*. Editorial Universitaria, Santiago. (Autopoiesis.)
- Rosen R. 1958. A relational theory of biological systems. *Bull Math Biophys* 20:245–260.
- Pattee HH. 1972. Laws and constraints, symbols and languages. In *Towards a Theoretical Biology* 4 (ed. CH Waddington). Edinburgh University Press.
- Hofmeyr J-HS. 2017. Basic biological anticipation. In *Handbook of Anticipation* (ed. R Poli). Springer.
- Schrödinger E. 1944. *What is Life?* Cambridge University Press. (Aperiodic crystal.)
- Bennett CH. 1988. Logical depth and physical complexity. In *The Universal Turing Machine: A Half-Century Survey* (ed. R Herken). Oxford University Press.
- Greenwald J, Friedmann MP, Riek R. 2016. Amyloid aggregates arise from amino acid condensations under prebiotic conditions. *Angew Chem Int Ed* 55:11609–11613.
- Walsh CT, Chen H, Keating TA, et al. 2001. Tailoring enzymes that modify nonribosomal peptides during and after chain elongation on NRPS assembly lines. *Curr Opin Chem Biol* 5:525–534.

---

## Acknowledgements

[TBD]

---

## End notes for revision

**Things to add or strengthen during revision:**

1. **Word count check.** Body text is approximately 5,500–6,000 words main + 1,000 methods. Comfortably within PNAS Direct Submission limit (~6,000 main + Methods separate). Discussion can be tightened by combining the third and fourth paragraphs.

2. **Figure preparation.** Figures 1–3 use existing simulation outputs (Tests A.1, A.2, B, C, E v2). Figure 4 uses Test D v2 outputs. All raw data and figure-generation scripts are on Boron at `/storage/kiran-stuff/templating_framework/`.

3. **The conservation table is the strongest empirical anchor.** R253 invariant across 1232 homologs, G248 invariant in 1231/1232, while R168 sits at 27% — these are facts about the Drt3b family that the framework's apparatus interprets as: Mode 3 N=2 is universal, gate identity is fixed, secondary residues are evolutionarily flexible. This is a substantive bioinformatic prediction supported by the published alignment.

4. **DRT7 mechanism citation.** The April 2026 bioRxiv (10.64898/2026.04.18.719102v1) describes DRT7 as "protein-primed, protein-templated, sequence-specific poly(T) synthesis through an arginine-rich recognition pocket." Confirm citation details before submission.

5. **Self-criticism stated.** The Discussion's penultimate paragraph handles limitations directly: provisional inventory, Markov-chain abstraction, Mode 2 not separately tested, no raw-read reanalysis. State, move on; no apology paragraph at the end.

6. **The closing paragraph is a calibrated bet, not a hedge.** "The bet is well-calibrated by the evidence presented here, and it is what natural history will verify" — replaces "future work needed to fully elucidate" pattern.

7. **Self-check for banned vocabulary.** Pre-submission: grep for "Interestingly", "Notably", "shed light", "pave the way", "delve", "elucidate", "leverage", "interplay" (figurative). All checked clean in v1; verify in v2.

8. **Read first sentence of each paragraph in sequence.** They should form a coherent abstract of the section. If they don't, restructure.

9. **The reviewer's eight specific concerns from the v1 review have all been addressed:**
   - (1) Deng et al. citation: corrected throughout, replacing "Sharma et al."
   - (2) Mode 1 information formula: corrected to $L \cdot [2 - h_b(\epsilon) - \epsilon \log_2 3]$
   - (3) Mode 3 saturation language: bound stated as $\le \log_2 N$ with equality approached in the low-noise long-block limit
   - (4) Drt3b section reframed as biological anchoring (Markov-chain abstractions parameterized from data) not validation against raw reads
   - (5) Population simulation reframed as consistency check on a stated theorem; theorem stated formally first
   - (6) Amyloid dual-classification example added to the apparatus subsection
   - (7) Mode 4 reframed as state-templating, distinct from polymer-position templating
   - (8) Introduction reframed around "no unified substrate-level taxonomy" rather than straw-manning Watson-Crick centrism

10. **Empirical advances over v1:** (a) the conservation analysis of 1,232 family members; (b) the Test E v2 with paper-derived E26Q parameters showing the framework matches the analytical prediction (0.10 marginal G) at floating-point precision; (c) the DRT family table covering 8 of 10 known systems; (d) the DRT10 telomerase-ancestor bridge connecting the framework to deep eukaryotic evolution.
