| | |
|:---|:---|
| **Journal** | BioSystems (Elsevier) |
| **ISSN** | 0303-2647 |
| **CiteScore** | 4.0 |
| **Impact Factor** | 1.9 |
| **Article type** | Perspectives |
| **Open Access APC** | USD 3,050 (optional) |
| **Submission portal** | [https://www.editorialmanager.com/biosystems](https://www.editorialmanager.com/biosystems) |
| **Guide for Authors** | [https://www.sciencedirect.com/science/journal/03032647/publish/guide-for-authors](https://www.sciencedirect.com/science/journal/03032647/publish/guide-for-authors) |

---

# The Prebiotic Roots of Alpha-Helical Dominance: A Structural Bioinformatics Perspective on Protein Origins

**Authors:** Madson A. de Luna-Aragao<sup>1</sup>, Joao P. Bezerra Neto<sup>4</sup>, Denys E. da Silva Santos<sup>3</sup>, Ana M. Benko-Iseppon<sup>2</sup>

**Affiliations:**
- <sup>1</sup> Instituto de Ciencias Biologicas, Universidade Federal de Minas Gerais, Belo Horizonte, Brazil
- <sup>2</sup> Departamento de Genetica, Universidade Federal de Pernambuco, Recife, Brazil
- <sup>3</sup> Departamento de Quimica Fundamental, Universidade Federal de Pernambuco, Recife, Brazil
- <sup>4</sup> Instituto de Ciencias Biologicas, Universidade de Pernambuco, Recife, Brazil

**Corresponding authors:** Madson A. de Luna-Aragao, madsondeluna@ufmg.br; Ana M. Benko-Iseppon, ana.iseppon@ufpe.br

---

## Highlights

- Alpha-helices may carry a compositional record of prebiotic chemistry detectable across all domains of life
- Structural bioinformatics now enables systematic, cross-lineage study of protein evolutionary origins
- Ancestral sequence reconstruction, AlphaFold, and synthetic prebiotic peptide libraries open new avenues
- We identify five open questions and outline a convergent research agenda for the field
- Helices in AMPs, TFs, and ion channels link early structural chemistry to essential life functions

---

## Abstract

The alpha-helix is the dominant secondary structure element across all known proteomes and
is the structural basis of protein families essential for life maintenance, including
antimicrobial peptides, transcription factors, and ion channels. Whether this dominance
reflects a contingent outcome of Darwinian selection or a physicochemical predisposition
rooted in the prebiotic amino acid inventory has remained largely unexplored at
quantitative scale. Large-scale structural bioinformatics of the CATH S40 non-redundant
mainly-alpha dataset (7,997 domains, 117,665 helical segments, 2,368,790 residues) reveals
that the two strongest helix-forming residues (Alanine and Leucine, both canonical
prebiotic amino acids) jointly occupy approximately 24.0% of all alpha-helical positions,
and that this enrichment pattern is reproduced across 4 domains of life, 17 kingdoms,
74 phyla, and 1,655 unique taxa. These observations motivate a broader research agenda
at the intersection of structural biology, prebiotic chemistry, and evolutionary
bioinformatics. Here we present our view of the open questions in this field, the
conceptual and technical challenges that have limited progress, and the emerging
approaches (ancestral sequence reconstruction, generative structural AI, synthetic
prebiotic peptide chemistry, and membrane-interface biophysics) that are now positioned
to address them. We argue that the alpha-helix represents an exceptional model system
for understanding how physicochemical constraints and evolutionary forces jointly shape
the molecular record of life, and we outline a convergent research agenda for the next
decade.

**Keywords:** alpha-helix; prebiotic chemistry; origins of life; protein evolution;
ancestral sequence reconstruction; structural bioinformatics; molecular evolution

---

## 1. Setting the Scene: The Alpha-Helix as a Window into Early Molecular Evolution

Among all the recurring patterns in the molecular world of life, few are as universal
and as structurally precise as the alpha-helix. Described by Pauling and Corey in 1951
[1], the alpha-helix is stabilized entirely by local backbone hydrogen bonds between the
carbonyl oxygen of residue *i* and the amide nitrogen of residue *i*+4. This
self-sufficient, intramolecular stabilization mechanism means that helical stability is
determined primarily by the amino acid sequence of the helix itself, without requiring
the long-range tertiary contacts that govern beta-sheet formation and folded globular
structure [2]. A peptide of as few as 8-15 residues can adopt a stable helical
conformation in solution when composed of appropriate amino acids [3]. This thermodynamic
accessibility sets the alpha-helix apart from all other secondary structure elements as
the most plausible candidate for the earliest functional polypeptide structure.

The biological universality of the alpha-helix is equally striking. Across all known
proteomes, helical residues account for approximately 32% of all positions in solved
structures. The "all-alpha" and "alpha/beta" CATH classes together represent the majority
of known protein domain folds [4,5]. More significantly, a disproportionate number of
protein families with irreplaceable functions in cellular life are predominantly or
exclusively helical. Antimicrobial peptides (AMPs), front-line innate immune effectors
conserved across all domains of life, adopt amphipathic alpha-helical conformations as
their functional state, inserting into bacterial membranes through a mechanism that
requires minimal sequence complexity [6,7]. Transcription factors of the helix-turn-helix,
basic helix-loop-helix, and leucine zipper families are helical at their DNA-binding and
dimerization interfaces [8]. Voltage-gated ion channels, Na+/K+-ATPases, and
mechanosensitive channels are built primarily from transmembrane alpha-helices [9]. The
deep convergence of essential cellular functions on the alpha-helical fold is not a
coincidence to be explained away: it is a biological signal demanding mechanistic
interpretation.

Yet the question of why the alpha-helix occupies this privileged position has received
surprisingly little systematic attention from a deep evolutionary perspective. The
standard answer, that evolution repeatedly discovered and retained the helix because it
is functionally useful, is circular: it explains the helix's prevalence by its
prevalence. A more informative question is whether there is a physicochemical reason,
prior to or independent of the specific functional contexts of modern proteins, why the
helix is so abundant and so compositionally constrained. We believe the answer to this
question lies in the intersection of prebiotic chemistry, protein thermodynamics, and
deep evolutionary bioinformatics, and that this intersection is now accessible to
systematic investigation.

---

## 2. What the Data Are Already Telling Us

Before identifying the open questions, it is worth establishing what is already
measurable. Structural bioinformatics at the scale of the CATH S40 non-redundant
mainly-alpha dataset provides an initial quantitative anchor.

### 2.1 Prebiotic amino acids dominate alpha-helical composition

![Figure 1](imgs/helix_propensities.png)

**Figure 1. Amino acid helical propensity and prebiotic classification.** Log2-transformed
observed helical propensity (blue bars) versus Chou-Fasman theoretical propensity (orange
bars) for all 20 standard amino acids, computed from 2,368,790 residues across 7,997
non-redundant mainly-alpha domains. Dashed line: neutral propensity (P = 1.0). Prebiotic
amino acids are indicated by asterisks.

Across 2,368,790 residues from 7,997 non-redundant mainly-alpha protein domains, Alanine
(observed propensity 1.300; proteome enrichment relative to Swiss-Prot human reference:
1.550) and Leucine (propensity 1.255; enrichment 1.365) are the two strongest
helix-forming residues. Both are canonical prebiotic amino acids, consistently recovered
in Miller-Urey spark discharge experiments [10], meteoritic analyses of the Murchison
and Murray chondrites [11,12], and hydrothermal vent synthesis studies [13]. Together
they occupy approximately 24.0% of all alpha-helical residue positions across 90,263
alpha-helical segments. At the other end of the spectrum, Glycine (enrichment 0.474) and
Proline (enrichment 0.379), also prebiotically abundant, are the most depleted residues
in helical regions, excluded by backbone dihedral and hydrogen bond geometry
respectively [14,15].

**Table 1. Proteome enrichment of key amino acids in helical regions and their prebiotic status.**

| AA | Prebiotic? | Enrichment vs. human proteome | Structural basis |
|----|:----------:|------------------------------:|:----------------|
| A  | Yes | **1.550** | Small methyl side chain favors helical phi/psi |
| L  | Yes | **1.365** | Isobutyl side chain packs helical cores |
| E  | Yes | 1.286 | Charged cap and stabilizing interactions |
| I  | Yes | 1.213 | Branched aliphatic, hydrophobic core |
| G  | Yes | 0.474 | No side chain; wide phi/psi space, helix-breaking |
| P  | Yes | 0.379 | Cyclic, backbone NH absent; kink-introducing |

This bifurcation within the prebiotic amino acid set, between helix-enriched (A, L, I,
V) and helix-depleted (G, P) residues, is precisely predicted by backbone geometry and
is reproduced across 4 domains of life, 17 kingdoms, 74 phyla, and 1,655 unique taxa
(Figure 2). The pattern therefore cannot be an artifact of lineage-specific evolutionary
history.

### 2.2 Positional preferences reflect physical chemistry

![Figure 2](imgs/ncap_ccap_preferences.png)

**Figure 2. N-cap and C-cap residue preferences across 117,665 helical termini.**
Z-score normalized amino acid enrichment at positions N1-N3 (left) and C1-C3 (right)
of alpha-helices. The Asparagine enrichment at N2-N3 (z approximately +3.07) is fully
predicted by hydrogen bond geometry at unsatisfied backbone amides.

![Figure 3](imgs/shannon_entropy_heptad.png)

**Figure 3. Shannon entropy and Leucine frequency at each heptad position.** All seven
heptad positions reach 96.3-96.4% of the theoretical maximum entropy of 4.322 bits,
while Leucine maintains a consistent 10-11% frequency at every position regardless of
structural role.

Residue preferences at helix termini (N-cap Asparagine enrichment at z approximately
+3.07), heptad hydrophobicity (~39% hydrophobic fraction at all seven positions
independently of structural role), and the cross-positional Leucine excess (~10-11% at
every heptad position, Figure 3) are all consistent with physicochemical principles that
would operate in any peptide regardless of its evolutionary history. Near-maximum Shannon
entropy (96.3-96.4% of the 4.322-bit theoretical maximum; Figure 3) demonstrates that
evolution has broadly diversified composition at every heptad position, yet has not
erased the Leucine baseline signal.

These observations are informative and reproducible. What they cannot tell us, by
themselves, is causality, directionality, or mechanism. That is the terrain of the open
questions ahead.

---

## 3. Open Questions

We identify five interconnected open questions that we believe define the frontier of
this research area and that, together, constitute the most important conceptual
challenges for the next decade.

### OQ1: Is the prebiotic signal in modern helices a direct inheritance or a convergent outcome?

The cross-lineage universality of the Ala/Leu enrichment pattern (Figure 2, Section 2.1)
is consistent with two fundamentally different historical scenarios. In the **inheritance
scenario**, a compositional bias toward Ala and Leu was established in the earliest
helical peptides of the last universal common ancestor (LUCA) or its predecessors and
has been vertically transmitted, with some erosion, to all extant lineages. In the
**convergence scenario**, independent lineages, operating under the same backbone
geometry and hydrophobic exclusion constraints, independently enriched their helical
regions in Ala and Leu through parallel evolutionary trajectories.

These two scenarios make different predictions about ancestral helical composition and
about the correlation between phylogenetic depth and compositional divergence. Ancestral
sequence reconstruction approaches [16] applied to deeply conserved mostly-helical
protein families (ferredoxins, cytochrome subunits, ATPase subunits) could, in
principle, distinguish between them: under the inheritance scenario, reconstructed
ancestral sequences should show monotonically increasing Ala+Leu content toward the
root of the tree of life, while under convergence, no such systematic trend is predicted.
This question has profound implications for our understanding of LUCA's molecular
composition and for the debate over whether a protein world or an RNA world preceded
the common ancestor [17].

### OQ2: What was the role of membrane interfaces in selecting early helical composition?

The distribution of Eisenberg hydrophobic moments across 90,263 alpha-helical segments
(Figure 4) reveals a substantial population of strongly amphipathic helices (muH > 0.40),
consistent with membrane-interfacial or helix-bundle packing contexts. Membrane-first
models of the origin of life [18,19] propose that the earliest functional peptides were
membrane-active, analogous to modern AMPs. This is a compelling scenario: AMP function
requires minimal sequence complexity (an amphipathic helix that segregates hydrophobic
and charged residues) and is exactly the type of function that could emerge spontaneously
from a prebiotic peptide mixture composed of Ala, Leu, Glu, and Lys [6,7].

![Figure 4](imgs/hydrophobic_moment_distribution.png)

**Figure 4. Distribution of Eisenberg hydrophobic moments across 90,263 alpha-helical
segments.** Broad distribution with a notable high-muH tail (muH > 0.40), consistent
with a substantial population of amphipathic helices suited for membrane-interfacial
contexts.

If the earliest helical peptides were membrane-active, then the physical chemistry of
the lipid-water interface, which strongly favors Ala and Leu at the hydrophobic face
and Lys, Glu at the hydrophilic face, would have been the first selection pressure on
helical composition, prior to and independent of Darwinian sequence optimization. This
is, in our view, one of the most tractable open questions in the field: it makes specific
experimental predictions about prebiotic peptide behavior at lipid membranes that are
accessible to modern biophysical methods (see Section 5).

### OQ3: How did the genetic code reinforce or reshape the prebiotic helical signature?

Alanine (4 synonymous codons) and Leucine (6 synonymous codons) are among the most
degenerate amino acids in the standard genetic code. If the genetic code was established
after the prebiotic helical bias was already in place, codon degeneracy may have served
to reinforce and maintain the compositional signal in translated proteomes. If, on the
other hand, the genetic code was established in a context where helical peptides were
already functionally important, structural constraints may have influenced codon
assignment, as a co-evolutionary feedback between chemistry, structure, and coding.

![Figure 5](imgs/codon_degeneracy_vs_propensity.png)

**Figure 5. Codon degeneracy versus observed helical propensity for all 20 standard
amino acids.** Alanine and Leucine occupy the high-degeneracy, high-propensity quadrant.
Glycine and Proline, also with 4 synonymous codons, fall at the bottom of the propensity
axis, demonstrating that degeneracy alone does not explain helical enrichment.

This connects directly to the long-standing debate over the origin and evolution of the
genetic code [20,21]. The observation that codon degeneracy does not explain propensity
differences within the prebiotic set (Figure 5): Gly and Pro have 4 codons each yet are
strong helix-breakers, which suggests that the code does not create the bias but may
perpetuate it. Whether this reflects a causal relationship or an independent coincidence
remains open. Computational co-evolutionary models of the genetic code that incorporate
structural constraints [22] represent one path toward addressing this question.

### OQ4: What distinguishes the helix as a structural attractor from other secondary structures?

![Figure 6](imgs/umap_aa_composition.png)

**Figure 6. UMAP projection of the amino acid composition space of 117,665 helical
segments.** Continuous compositional manifold with Leucine- and Alanine-dominated
helices forming the highest-density central region, consistent with a compositional
attractor of the alpha-helical fold space.

The UMAP projection of per-helix amino acid composition (Figure 6) reveals a continuous
compositional manifold rather than discrete clusters, with Ala- and Leu-dominated
helices forming the dense central attractor. What makes this attractor specifically
associated with the alpha-helix, and not with beta-strands or other secondary structures?
Equivalent analyses applied to beta-strand-enriched protein classes would reveal whether
the prebiotic bias is helix-specific or whether each secondary structure type has its own
prebiotic attractor consistent with its backbone geometry. The alpha-helix and beta-strand
impose different steric and hydrogen bonding constraints on their constituent amino acids:
we predict that beta-strand regions should show enrichment in Val and Ile (also prebiotic)
and depletion of Ala, while the alpha-helical pattern shows the reverse. If confirmed,
this would indicate that the prebiotic amino acid set was not only abundant but was
differentiated in ways that mapped onto the full spectrum of secondary structure
propensities. The structural universe of proteins may, in this sense, be a direct
readout of prebiotic amino acid chemistry.

### OQ5: Can we reconstruct the compositional trajectory of helices from prebiotic peptides to LUCA?

The deepest open question is whether it is possible, in principle and in practice, to
reconstruct the compositional history of the alpha-helix from its prebiotic origins to
the modern protein universe. This would require integrating: (i) experimental data on
prebiotic peptide synthesis and helix formation; (ii) ancestral sequence reconstruction
of ancient helical protein families; (iii) comparative compositional analysis across all
domains of life; and (iv) phylogenetic models of amino acid substitution that explicitly
account for backbone structural constraints. No such integrated framework currently
exists, but the emergence of high-quality, large-scale structural databases (CATH,
AlphaFold Database [23]), improved ancestral reconstruction methods [16], and the
explosion of sequenced genomes across the tree of life makes this goal more achievable
than at any previous time.

---

## 4. Conceptual and Technical Challenges

Before the research agenda outlined in Section 5 can be fully executed, several
conceptual and technical challenges must be acknowledged.

### 4.1 Confounding effects of functional selection

The most serious confounding factor in any attempt to identify a prebiotic signal in
modern protein composition is that billions of years of functional selection have
substantially shaped helical composition independently of any prebiotic baseline.
Membrane-embedded transmembrane helices are enriched in hydrophobic residues by
functional necessity [9]; coiled-coil leucine zippers are enriched in Leu at heptad
position d by structural requirement [24]; AMPs are enriched in Lys and Arg by the need
to interact with anionic bacterial membranes [7]. Disentangling the prebiotic baseline
from functional signals requires careful stratification by protein family, subcellular
localization, and functional context. The near-maximum Shannon entropy at all heptad
positions (96.3-96.4% of the theoretical maximum) suggests that functional selection has
broadly diversified composition without eliminating the prebiotic floor, but a systematic
decomposition of the observed enrichment into prebiotic and functional components has not
yet been achieved.

### 4.2 Uncertainty in the prebiotic amino acid inventory

The prebiotic classification of amino acids used in structural bioinformatics studies is
necessarily based on experimental consensus from Miller-Urey experiments [10], meteoritic
analyses [11,12], and hydrothermal synthesis [13], but the actual abundance ratios in
prebiotic environments were almost certainly geochemically heterogeneous and temporally
variable. The abundance of Glycine relative to Alanine, for example, varies substantially
depending on the energy source, the pH, and the specific mineral environment [25,26].
This uncertainty limits the precision with which prebiotic predictions can be made and
means that the "prebiotic amino acid set" is itself an approximation that requires
explicit treatment of compositional uncertainty in future analyses.

### 4.3 The causality problem

Perhaps the deepest conceptual challenge is that the observed correlation between
prebiotic amino acid abundance and helical enrichment does not, by itself, establish
causality. There are at least three non-exclusive explanations: (i) prebiotic chemistry
established the helical bias, which was inherited by all subsequent life; (ii)
evolutionary convergence under the same physical constraints independently produced the
same pattern in each lineage; or (iii) the genetic code independently reinforced the
pattern by assigning more codons to residues with high helical propensity, for reasons
unrelated to early peptide chemistry. These explanations require different experimental
and computational approaches to distinguish and may not be mutually exclusive.

### 4.4 The pre-LUCA gap

The most significant challenge is perhaps that there is, by definition, no sequence or
structural record from the period between the first prebiotic peptides and LUCA. The
only access to this period is indirect: through the amino acid compositions of proteins
that LUCA must have possessed, inferred from the universal core of conserved proteins
[27,28], and through experimental reconstruction of prebiotic peptide chemistry. Bridging
this gap, connecting the chemistry of the prebiotic world to the molecular record of
LUCA, remains the central challenge of the field.

---

## 5. Emerging Approaches and Future Research Directions

Despite these challenges, several convergent technological and conceptual developments
position the field to make substantial progress in the coming decade. We outline five
research directions that we consider most promising.

### 5.1 Ancestral sequence reconstruction at deep phylogenetic scales

Maximum-likelihood and Bayesian ancestral sequence reconstruction [16,29] applied to
the universal core of deeply conserved mostly-helical proteins (including ferredoxins,
cytochromes, ATPase subunits, and signal recognition particle components) offers a
direct window into helical composition at or near LUCA. If reconstructed ancestral
sequences show monotonically increasing Ala+Leu content toward the root of the tree of
life, this would be the most direct evidence that the prebiotic helical bias was present
at LUCA and has been gradually eroded by evolutionary diversification. If no such trend
is observed, the convergence or genetic code hypotheses would be supported. This approach
is now technically feasible with available computational resources and the expanding
availability of deep-branching genomes from metagenomics [30,31].

### 5.2 Generative structural AI and compositional exploration

The AlphaFold revolution [23] and the emergence of protein language models (ESM-2 [32],
ProtTrans [33]) trained on hundreds of millions of protein sequences open a new
analytical space for this question. These models implicitly learn the statistical
regularities of protein sequences, including compositional biases by secondary structure
type. Probing these models for the degree to which they have internalized a prebiotic
compositional prior, using techniques from mechanistic interpretability and concept
activation vectors [34], could reveal whether the prebiotic helical bias is a structural
regularity learned from the data rather than an artifact of framing. More ambitiously,
generative protein design models could be used to explore the sequence space of
prebiotic-like helical peptides: generating sequences constrained to the prebiotic amino
acid inventory and evaluating their folding propensity using structure prediction, without
any reliance on evolutionary homology. This would allow the prebiotic question to be
addressed in silico at scale.

### 5.3 Experimental prebiotic peptide libraries and biophysical characterization

The ultimate experimental approach to the prebiotic helical question is to synthesize
peptide libraries from a prebiotically realistic amino acid mixture and measure their
structural properties. Solid-phase peptide synthesis now allows the routine generation
of random peptide libraries of defined length and composition [35]. A library of 15-25
residue peptides drawn from an amino acid mixture weighted by prebiotic abundance
estimates could be characterized by circular dichroism (CD) spectroscopy for helical
content in aqueous solution and in membrane-mimetic environments (TFE, DPC micelles,
lipid vesicles). If prebiotic-weighted libraries show significantly higher helical
content than equimolar controls, this would provide direct experimental support for the
idea that early chemistry was biased toward helical structures. Library members with
highest helical content could subsequently be tested for membrane activity, providing a
direct link between compositional chemistry and functional potential.

### 5.4 Transmembrane and amphipathic helix analysis as a model subsystem

Transmembrane (TM) helices represent a distinct subsystem within the broader alpha-helical
universe: they are defined by their membrane-embedded context, they interact predominantly
with the lipid bilayer rather than with other proteins, and they are functionally
essential in all domains of life. If the prebiotic bias is driven in part by
membrane-interface physics, TM helices should show a stronger Ala/Leu enrichment than
soluble helices of equivalent length. Systematic comparison of TM and soluble helical
composition using OPM- or PDBTM-annotated datasets [36], stratified by functional class
and phylogenetic depth, would test this prediction directly and would connect the
structural bioinformatics approach to the membrane-first origin-of-life framework.

### 5.5 Integration with RNA-world and RNA-protein coevolution models

The relationship between the prebiotic helical bias and the RNA world remains largely
unexplored. If the earliest peptides were synthesized on an RNA template (or on a
proto-ribosome), the compositional bias of those peptides would have been determined
partly by the coding capacity of the earliest RNA sequences and partly by the
selectivity of aminoacyl-tRNA synthetase precursors. Models of RNA-protein coevolution
[37,38] that incorporate structural constraints on the encoded peptides, specifically the
preference for helical conformations, could reveal whether RNA sequences were
pre-adapted to encode helical peptides or whether the structural bias emerged through
subsequent optimization. The observation that Alanine and Leucine also carry among the
highest codon degeneracies in the modern genetic code (Figure 5) raises the question of
whether this degeneracy is causally connected to the structural role of these residues in
early peptide evolution.

---

## 6. A Convergent Research Agenda

The five research directions outlined in Section 5 are not independent: their power lies
in their convergence. We envision a research agenda in which structural bioinformatics
provides the quantitative baseline (what is the compositional pattern in modern
proteins?), ancestral reconstruction traces the pattern back toward LUCA (when and how
did it arise?), experimental prebiotic peptide chemistry tests functional plausibility
(could early peptides with this composition be functional?), and generative AI explores
the full sequence space consistent with prebiotic constraints (how broad is the prebiotic
helical landscape?).

![Figure 7](imgs/taxonomy_cladogram.png)

**Figure 7. Taxonomic breadth of the structural dataset.** Circular cladogram spanning
4 domains of life, 17 kingdoms, 74 phyla, and 1,655 unique taxa. The reproducibility
of the prebiotic helical signal across this breadth (Tables 2-3) establishes the
empirical baseline for the comparative and ancestral reconstruction analyses proposed in
Section 5.

**Table 2. Taxonomic coverage of the structural dataset supporting the cross-lineage claim.**

| Domain of life | n(PDBs) | n(taxa) | n(kingdoms) | n(phyla) |
|---|---:|---:|---:|---:|
| Eukaryota | 4,061 | 406 | 3 | 30 |
| Bacteria | 3,720 | 894 | 4 | 25 |
| Viruses | 373 | 208 | 6 | 15 |
| Archaea | 424 | 84 | 4 | 7 |
| **Total** | **8,265** | **1,655** | **17** | **74** |

**Table 3. Core compositional observations motivating the research agenda.**

| Observation | Value | Scope | Structural basis |
|---|---|---|---|
| Ala enrichment in alpha-helices vs. human proteome | 1.550 | 90,263 segments | Methyl side chain favors helical phi/psi |
| Leu enrichment in alpha-helices vs. human proteome | 1.365 | 90,263 segments | Hydrophobic core packing |
| Ala+Leu fraction of all alpha-helical positions | ~24.0% | 90,263 segments | Combined effect |
| Gly depletion in alpha-helices | 0.474 | 90,263 segments | Wide phi/psi; helix-incompatible |
| Pro depletion in alpha-helices | 0.379 | 90,263 segments | Cyclic; no backbone NH |
| Asn enrichment at N2 cap (z-score) | +3.07 | 117,665 termini | H-bond to unsatisfied backbone NH |
| Heptad Shannon entropy (mean) | 4.16 bits (96.3% max) | 7 positions | Near-uniform composition |
| Leu frequency at all heptad positions | ~10-11% | 7 positions | Position-invariant prebiotic floor |
| Taxa spanning the pattern | 1,655 unique (4 domains) | 8,265 structures | Cross-lineage generality |

The most transformative contribution of this agenda would be a mechanistic, time-resolved
account of how the alpha-helical fold acquired its compositional identity: from the
prebiotic soup, through the RNA world, through LUCA, and into the remarkable structural
diversity of modern proteomes. This is not a question with a simple answer, but it is
now, for the first time, a question with tractable experimental and computational
approaches.

---

## 7. Implications for the Origins of Life Field

The perspective we have outlined has implications that extend beyond structural biology.
If the alpha-helix is indeed a structural fossil of prebiotic chemistry, retaining a
measurable record of the amino acid inventory of the early Earth, then the protein
structural databases are, in a real sense, a geological record of prebiotic chemistry,
written in the language of amino acid frequencies rather than in the language of isotope
ratios or mineral assemblages. Reading this record requires the same rigor, and the same
caution about alternative interpretations, that is applied to any other geological archive.

This has practical implications for how we think about the origin of life. The question
"how did the first proteins arise?" has traditionally been addressed either through
experimental prebiotic chemistry (what peptides can form abiotically?) or through
phylogenetics (what proteins did LUCA possess?). The structural bioinformatics approach
we advocate adds a third line of evidence: what does the composition of modern protein
structures tell us about the physicochemical constraints that shaped the earliest peptides?
All three lines of evidence are imperfect and incomplete, but their convergence or
divergence will be informative.

We also note that the alpha-helix's role in essential proteins (AMPs, ion channels,
transcription factors) suggests that helical architecture may have been functionally
central from very early in the history of life, not merely a structural scaffold but an
active participant in the emergence of membrane function, gene regulation, and molecular
recognition. Understanding the compositional evolution of the helix is therefore not only
a question about structural biology: it is a question about the chemical origins of
biological function itself.

---

## 8. Conclusions

The alpha-helix stands at the crossroads of prebiotic chemistry and modern molecular
evolution. Structural bioinformatics at the scale of thousands of non-redundant protein
domains reveals a reproducible, cross-lineage compositional pattern in which prebiotic
amino acids with high helical propensity are systematically enriched, and those with low
propensity are systematically depleted, in a manner predicted by backbone geometry and
consistent with physicochemical constraints that would have operated in prebiotic
peptides. The pattern is visible across 4 domains of life, 17 kingdoms, and 74 phyla,
suggesting it is not a lineage-specific artifact but a deep signature of protein
structural history.

Yet the most important questions remain open. Is this signal a direct inheritance from
prebiotic chemistry or a product of convergent evolution under shared physical
constraints? What role did membrane interfaces play in selecting early helical
composition? How did the genetic code reinforce or reshape this prebiotic bias? Can we
reconstruct the compositional trajectory of the helix from the early Earth to LUCA and
beyond?

We argue that these questions are now tractable, and that their answers require a
convergent research agenda combining ancestral sequence reconstruction, generative
structural AI, experimental prebiotic peptide chemistry, and transmembrane helix
biophysics. The alpha-helix, as the most thermodynamically accessible, functionally
essential, and compositionally distinctive secondary structure element, represents an
exceptional model system for addressing one of the deepest questions in biology: what
chemical logic shaped the molecular architecture of the first living systems?

---

## Data Availability

Structural data: CATH database (v4.3; https://www.cathdb.info).
PDB files: RCSB Protein Data Bank (https://www.rcsb.org).
Taxonomy data: NCBI Entrez API (https://www.ncbi.nlm.nih.gov/entrez).
Reference proteome: UniProt Swiss-Prot human canonical sequences (https://www.uniprot.org).
Analysis code and processed dataset available at [repository URL] upon acceptance.

---

## Author Contributions

[CRediT statement to be completed]

---

## Declaration of Competing Interests

The authors declare no competing interests.

---

## Funding

[Funding statement to be completed]

---

## Acknowledgements

[To be completed]

---

## References

[1] Pauling L, Corey RB, Branson HR. The structure of proteins: two hydrogen-bonded
helical configurations of the polypeptide chain. *PNAS*. 1951;37:205-211.

[2] Creighton TE. *Proteins: Structures and Molecular Properties*. 2nd ed. New York:
W.H. Freeman; 1993.

[3] Marqusee S, Robbins VH, Baldwin RL. Unusually stable helix formation in short
alanine-based peptides. *PNAS*. 1989;86:5286-5290.

[4] Sillitoe I, et al. CATH: increased structural coverage of functional space.
*Nucleic Acids Research*. 2021;49(D1):D266-D273.

[5] Andreeva A, et al. SCOP2 prototype: a new approach to protein structure mining.
*Nucleic Acids Research*. 2014;42(D1):D310-D314.

[6] Zasloff M. Antimicrobial peptides of multicellular organisms. *Nature*.
2002;415:389-395.

[7] Brogden KA. Antimicrobial peptides: pore formers or metabolic inhibitors in bacteria?
*Nature Reviews Microbiology*. 2005;3:238-250.

[8] Luscombe NM, Austin SE, Berman HM, Thornton JM. An overview of the structures of
protein-DNA complexes. *Genome Biology*. 2000;1:reviews001.

[9] Doyle DA, et al. The structure of the potassium channel: molecular basis of K+
conduction and selectivity. *Science*. 1998;280:69-77.

[10] Miller SL. A production of amino acids under possible primitive Earth conditions.
*Science*. 1953;117:528-529.

[11] Cronin JR, Pizzarello S. Enantiomeric excesses in meteoritic amino acids. *Science*.
1997;275:951-955.

[12] Botta O, Bada JL. Extraterrestrial organic compounds in meteorites. *Surveys in
Geophysics*. 2002;23:411-467.

[13] Huber C, Wachtershäuser G. Peptides by activation of amino acids with CO on (Ni,Fe)S
surfaces. *Science*. 1998;281:670-672.

[14] Ramachandran GN, Ramakrishnan C, Sasisekharan V. Stereochemistry of polypeptide
chain configurations. *Journal of Molecular Biology*. 1963;7:95-99.

[15] MacArthur MW, Thornton JM. Influence of proline residues on protein conformation.
*Journal of Molecular Biology*. 1991;218:397-412.

[16] Yang Z. PAML 4: phylogenetic analysis by maximum likelihood. *Molecular Biology and
Evolution*. 2007;24:1586-1591.

[17] Orgel LE. The origin of life: a review of facts and speculations. *Trends in
Biochemical Sciences*. 1998;23:491-495.

[18] Deamer D. The first living systems: a bioenergetic perspective. *Microbiology and
Molecular Biology Reviews*. 1997;61:239-261.

[19] Szostak JW, Bartel DP, Luisi PL. Synthesizing life. *Nature*. 2001;409:387-390.

[20] Koonin EV, Novozhilov AS. Origin and evolution of the genetic code: the universal
enigma. *IUBMB Life*. 2009;61:99-111.

[21] Trifonov EN. Consensus temporal order of amino acids and evolution of the triplet
code. *Gene*. 2000;261:139-151.

[22] Sengupta S, Higgs PG. Pathways of genetic code evolution in ancient and modern
organisms. *Journal of Molecular Evolution*. 2015;80:229-243.

[23] Jumper J, et al. Highly accurate protein structure prediction with AlphaFold.
*Nature*. 2021;596:583-589.

[24] Lupas AN, Bassler J, Dunin-Horkawicz S. The structure and topology of alpha-helical
coiled coils. *Subcellular Biochemistry*. 2017;82:95-129.

[25] Parker ET, et al. Primordial synthesis of amines and amino acids in a 1958 Miller
H2S-rich spark discharge experiment. *PNAS*. 2011;108:5526-5531.

[26] Kebukawa Y, et al. Organic matter in the Belgica-7904 carbonaceous chondrite.
*Meteoritics and Planetary Science*. 2010;45:394-408.

[27] Lecompte O, et al. Comparative analysis of ribosomal proteins in complete genomes:
an example of reductive evolution at the domain scale. *Nucleic Acids Research*.
2002;30:5382-5390.

[28] Fournier GP, Gogarten JP. Rooting the ribosomal tree of life. *Molecular Biology
and Evolution*. 2010;27:1792-1801.

[29] Hochberg GKA, Thornton JW. Reconstructing ancient proteins to understand the causes
of structure and function. *Annual Review of Biophysics*. 2017;46:247-269.

[30] Hug LA, et al. A new view of the tree of life. *Nature Microbiology*. 2016;1:16048.

[31] Parks DH, et al. A standardized bacterial taxonomy based on genome phylogeny
substantially revises the tree of life. *Nature Biotechnology*. 2018;36:996-1004.

[32] Lin Z, et al. Evolutionary-scale prediction of atomic-level protein structure with
a language model. *Science*. 2023;379:1123-1130.

[33] Elnaggar A, et al. ProtTrans: toward understanding the language of life through
self-supervised learning. *IEEE Transactions on Pattern Analysis and Machine
Intelligence*. 2022;44:7112-7127.

[34] Kim B, et al. Interpretability beyond classification: quantitative testing with
concept activation vectors (TCAV). *ICML Proceedings*. 2018.

[35] Merrifield RB. Solid phase peptide synthesis: I. The synthesis of a tetrapeptide.
*Journal of the American Chemical Society*. 1963;85:2149-2154.

[36] Lomize MA, et al. OPM database and PPM web server: resources for positioning of
proteins in membranes. *Nucleic Acids Research*. 2012;40(D1):D370-D376.

[37] Caetano-Anolles G, et al. The origin and evolution of the ribosome. *Biology
Direct*. 2011;6:43.

[38] Yarus M. The genetic code and RNA-amino acid affinities. *Life*. 2017;7:13.

---

## Submission Checklist (BioSystems Perspectives)

- [ ] Article type explicitly set to "Perspectives" at submission
- [ ] Cover letter explaining visionary/future-oriented scope and fit to BioSystems
- [ ] Title page: full author list, affiliations, corresponding author contact
- [ ] Abstract (~200-250 words), visionary tone
- [ ] Keywords (4-6 terms)
- [ ] Highlights (3-5 bullet points, max 85 characters each including spaces)
- [ ] Open questions clearly identified and framed as challenges for the field
- [ ] Future research directions section with actionable, specific proposals
- [ ] Figures at minimum 300 DPI for review; 600 DPI for line art in final submission
- [ ] Figure captions in manuscript body, not embedded in figures
- [ ] Tables as editable text, not images
- [ ] References in Elsevier numbered format [1], [2], ...
- [ ] Authorship declaration (all authors)
- [ ] Competing interests statement
- [ ] Generative AI use declared if applicable
- [ ] Funding sources disclosed
- [ ] Data availability statement
- [ ] No concurrent submission confirmed

---

*Working draft, Perspectives format for BioSystems (Elsevier).*
*Status: draft; all figures generated, text in progress.*
