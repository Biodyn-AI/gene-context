# Single-cell foundation models represent genes in context, but the information is co-expression

**Ihor Kendiukhov · 2026**

## Abstract

A single-cell foundation model reads a cell as a sequence of its genes and produces, at each layer, a contextual
representation of every gene — the gene's vector given the other genes in that cell. Whether that representation
carries usable information about a gene's role in a particular cell is unknown, and is the single-cell analogue
of word-sense disambiguation in language models. We give a systematic answer. Across three architectures
(MaxToki, scGPT, STATE) a gene's representation is reshaped by cellular context in a way that is gene-specific
(cross-cell agreement 0.759 versus 0.001 for a gene-shuffled null), reproducible across independent cells, not an
artefact of the gene's rank position (96% survives rank control), and replicable on an independent developmental
dataset. It is causally used: steering a functional direction into half a cell's genes shifts the model's own
predictions for the other genes (p = 7×10⁻⁴). It is learned rather than architectural: an untrained model of the
same architecture organises context far less (z = +3.8 versus +21 for the trained model). It also increases with
model size, though this scaling increment is partly confounded by embedding width. The context modulation is
organised along interpretable functional axes, but against size-matched co-expression nulls this organisation
does not robustly exceed co-expression — the strongest axis is only at the p = 0.05 border and the others are
non-significant — and two independent probes for novel context-dependent gene function are null. The
useful, ceiling-immune output is a model metric: contextualisation and its functional organisation differ
sharply by tokenisation (rank-based versus value-binned versus protein-sequence-token) and scale with size. We
report the phenomenon, place a boundary on what it encodes, and give the metric.

**Keywords:** single-cell foundation models; contextual gene representations; mechanistic interpretability;
co-expression; gene embeddings; word-sense disambiguation; model characterisation.

## 1. Introduction

A recurring theme in the interpretability of large language models is that a network trained on a narrow
predictive task builds internal representations of structure it was never taught. Sequence models of Othello and
chess build representations of the board and of latent variables such as player skill [@li2023emergent;
@nanda2023emergent; @karvonen2024emergent]; language models linearly encode space and time along identifiable
directions [@gurnee2024language]; and, in the best cases, these representations are not decorative but *causally
used* — a claim one tests by intervening on the representation and measuring the effect on behaviour
[@turner2023activation; @zou2023representation; @park2024linear].

A closely related line concerns *contextual* representations. Contextual word embeddings assign the same word a
different vector in every sentence, and the geometry of those vectors encodes word sense: "bank" near a river
and "bank" in finance occupy separable regions [@peters2018; @ethayarajh2019; @coenen2019]. A cottage industry
of methodological care surrounds these measurements, because representation spaces are strongly anisotropic (vectors cluster in a narrow cone rather than spreading
evenly, so any two look similar) and a handful of "rogue" dimensions can dominate cosine similarity [@ethayarajh2019; @muviswanath2018; @cai2021;
@timkey2021].

Single-cell foundation models (SCFMs) offer a clean biological version of the same question
[@geneformer2023; @scgpt2024; @scfoundation2024; @scbert2022]. They treat each gene as a token and each cell as
a sequence of its expressed genes, and are trained by masked or autoregressive prediction over
expression-ranked genes — in effect, learning which genes tend to be active together. Their gene *token*
embeddings (one fixed vector per gene) have been studied; their *contextual* gene representations (a gene's
vector inside a particular cell) are used as features downstream [@scfoundation2024; @bian2024] but have not
been characterised as objects in their own right. Yet this is exactly where a model could, in principle,
represent that a gene plays a different role in a macrophage than in a T cell — the biological analogue of word
sense, and, as W. Gilpin has noted, an open direction that the "cell as a token" framing invites but does not
test [@Gilpin2025CellToken].

We organise the question into two levels. **Level 1 (a model-fact):** is a gene's representation reshaped by
context in a structured, gene-specific way, and is that reshaping learned and used? **Level 2 (biology):** does
the reshaping reveal something *biologically novel* about the gene's role in each context — beyond what
co-expression already tells us? We answer Level 1 affirmatively and comprehensively, and Level 2 negatively and
with matched-null rigour. The negative is not a failure to report but the paper's second contribution: it
locates, precisely and with controls, the ceiling that a growing critical literature has found for SCFMs in
other guises [@kedzierska2025zeroshot; @boiarsky2024deeper; @ahlmanneltze2025perturbation; @kendiukhov2026attention]. The paper thus makes three contributions: (1) the Level-1 model-fact, characterised across architectures, scale, and a causal test; (2) the Level-2 co-expression ceiling, established with size-matched nulls; and (3) contextualisation as a model-characterisation metric that distinguishes tokenisation schemes (§4.7).

## 2. Related work

**Contextual representations and word sense.** Our methods descend directly from the NLP interpretability
literature on contextualised embeddings. Ethayarajh [@ethayarajh2019] showed that representations are anisotropic
in every layer, that same-word self-similarity falls in upper layers (they become more context-specific), and
that a static vector explains little of a word's contextual variance — and, crucially, that every cosine must be
corrected against an anisotropy baseline. Timkey and van Schijndel [@timkey2021] showed that a few rogue
dimensions dominate similarity and that a scalar correction is insufficient; Mu et al. [@muviswanath2018] and
Cai et al. [@cai2021] give the standard postprocessing and the local-isotropy picture. Coenen et al.
[@coenen2019] and Garí Soler and Apidianaki [@garisoler2021] show word senses form separable clusters and that
polysemy level is itself readable. We import both the measurements and the controls.

**Emergent, causally-used latent variables.** That a model represents and *uses* structure it was never given is
established for board games [@li2023emergent; @nanda2023emergent; @karvonen2024emergent] and for physical
coordinates in language models [@gurnee2024language]. We adopt this literature's evidentiary standard — a probe
shows only correlation, so a causal claim requires an intervention against a matched control — and its steering
machinery [@turner2023activation; @zou2023representation], of which the extraction and amplification of
interpretable feature directions in production models [@templeton2024scaling] is the large-scale instance. Our
causal result is a steer-here/read-there test in this mould.

**Single-cell foundation models and their representations.** Geneformer [@geneformer2023], scGPT [@scgpt2024],
scFoundation [@scfoundation2024], scBERT [@scbert2022], UCE [@uce2023] and STATE [@state2025] differ in
tokenisation (rank-ordered gene identity, value-binned expression, or protein-sequence-initialised gene tokens),
and the review of Bian et al. [@bian2024] defines the "contextual gene embedding" our work dissects. MaxToki
[@maxtoki2026], the primary model here, is a Llama-architecture, Geneformer-lineage model at 217M and 1B scales.
scGPT reports cell-state-specific gene networks via attention, the closest prior to a "context changes the gene"
claim — but operationalised as a pairwise network, never as a property of the representation vector, and without
a decomposition or a null. The *static* gene-embedding geometry of these models has been mapped — its spectral
structure [@kend2026spectral], the topological and geometric structure it does and does not learn
[@kend2026topo141], and a decomposition of what a protein language model already supplies versus what single-cell
pretraining adds [@kend2026twoaxes] — but the *contextual* representation, the object here, has not.

**What SCFMs do not capture.** A critical literature finds SCFMs often fail to beat simple baselines: on
zero-shot tasks [@kedzierska2025zeroshot], under deeper evaluation [@boiarsky2024deeper], on perturbation
prediction [@ahlmanneltze2025perturbation; @bendidi2024pca], and across utilities [@liu2026sceval]; and sparse
autoencoders find organised biological knowledge but minimal regulatory logic in these models
[@kend2026sae], while SCFM attention re-expresses co-expression rather than unique regulatory signal
[@kendiukhov2026attention]. Our
co-expression ceiling is the representation-geometry instance of this pattern, established with a matched null.

**The biology of the ceiling.** Co-expression is the canonical guilt-by-association signal
[@vanDam2018Coexpression], and gene-embedding methods from Gene2vec [@Du2019Gene2vec] to GenePT
[@Chen2023GenePT] are built on it. Genuine context-dependent gene function exists — moonlighting proteins
[@Jeffery1999Moonlighting], lineage-determining transcription factors that select different enhancers with
different partners by cell type [@Heinz2010LDTFs] — but much of it is post-transcriptional and so, we argue,
structurally invisible to an expression-only model, which is why our Level-2 tests are informative negatives.

## 3. Methods

**Extraction.** We run cells from Tabula Sapiens (immune, kidney, lung [@tabulasapiens2022]) through a model and
record, for each gene in a panel, its hidden state at that gene's token position, averaged over cells of one
cell type. The headline MaxToki-217M extraction covers 12 cell types × 1000 cells, a 6000-gene panel, taps at
layers 0, 4, 8, 11. Four design choices, each fixing a measured failure of a pilot: (1) **pairwise contexts**,
not the all-context intersection, which leaves ~65 housekeeping genes — the "stopwords" that are *most*
context-sensitive [@ethayarajh2019] — whereas pairwise comparison recovers ~5,700 genes per pair; (2)
**cell-level split halves**, each per-gene mean estimated twice from disjoint cells for an honest noise floor;
(3) an **occurrence cap** (50 tokens per gene, context, partition), because a per-gene mean has standard error
∝1/√n and unequal counts inject context-dependent noise into any interaction term; (4) **per-dimension
z-scoring** before any cosine, to remove the rogue-dimension inflation a scalar anisotropy correction misses
[@timkey2021].

**Decomposition.** For gene *g* in cell type *c*, write *v(g,c) = μ + a(g) + b(c) + E(g,c)*: gene main effect
*a(g)* (identity), context main effect *b(c)* (the **averaging null** — attention pools over the cell, so
everything drifts together), and gene-specific interaction *E(g,c)*, the only term where biology can live.

**Contextualisation strength (EXCESS).** For a cell-type pair (c₀,c₁), the *crowd-removed shift* of gene g is
δ(g) = [v(g,c₁) − v(g,c₀)] − mean_g′[v(g′,c₁) − v(g′,c₀)] — the gene's own shift minus the average shift over all
genes (the "crowd", i.e. the context main effect). We compute δ twice, from two disjoint halves of the cells
(δ₀, δ₁), and set EXCESS = mean_g cos(δ₀(g), δ₁(g)) − mean_g cos(δ₀(g), δ₁(perm g)): the same gene's shift
agreeing across independent cells, minus a gene-shuffled control. EXCESS = 0 means pure averaging (no
gene-specific context response); it is averaged over all cell-type pairs.

**Functional axes and nulls.** A functional axis u is the unit vector separating two gene classes (from gene
ontology) in the gene main effect: u ∝ mean(a(g)|class A) − mean(a(g)|class B). Its **reproducible interaction
power** is: project each gene's two-way interaction term E(g,c) onto u to get a scalar per (gene, context), then
average the product of the two cell-partitions' scalars over count-balanced (gene, context) cells — a covariance
that survives only if the modulation reproduces across independent cells. We compare this power for u against
(i) matched **random-axis** nulls (random gene groupings of the same pole sizes); (ii) **size- and
co-expression-matched** null gene-set pairs — the decisive control, since functional gene sets are co-expressed —
where co-expression coherence is the mean pairwise Pearson correlation of expression across the extraction cells
within the gene set, the null pairs are drawn to match the functional poles on both size and coherence, and the
functional axis's power is reported as an empirical percentile among them (skew-robust, because the null power
distribution is right-skewed); and (iii) a **representation-tightness-matched** null (tightness = mean cosine
among the axis genes' own representations). Full parameters in §7.

## 4. Results

Sections 4.1–4.3, 4.5 and 4.7–4.8 establish **Level 1** (the model-fact: contextualisation is real, not
position, replicable, causal, learned, and scaling); §4.4 and §4.6 test **Level 2** (novel biology) and find the
co-expression ceiling.

### 4.1 Genes are contextualised gene-specifically, and it peaks in the early-middle layers

The gene-specific context response is large and highly reproducible (Table 1, Fig. 1). At layer 4 the gene-specific EXCESS is **+0.758**: the same gene's crowd-removed shift agrees at cosine +0.759
across independent cells, while a *different* gene's
agrees at **+0.001** — the effect is entirely gene-specific. The context main effect replicates at +0.996,
confirming the measurement works. Layer 0 (the context-free embedding output) returns **exactly 0.000**, the
sanity check that the pipeline reports nothing when there is nothing. A scan across layers 0, 4, 8 and 11 places the peak at layer 4,
rising through the early layers and decaying toward the output.

| Layer | EXCESS (same − diff) | same-gene | diff-gene | main-effect replication |
|---|---|---|---|---|
| 0 (embedding) | **+0.000** | +0.000 | +0.000 | — |
| 4 | **+0.758** | +0.759 | +0.001 | +0.996 |
| 8 | +0.665 | +0.666 | +0.001 | +0.996 |
| 11 | +0.627 | +0.628 | +0.001 | +0.997 |

*Table 1. Gene-specific context response by layer (MaxToki-217M). 66 cell-type pairs; ~5,700 genes are shared
per pair, of which a median ~1,700 remain after the occurrence-cap count-balancing (§3); bootstrap CIs ±0.002.*

![**Figure 1.** Gene-specific context response by layer. Same-gene agreement across independent cells (blue) far exceeds the gene-shuffled null (grey); layer 0, the context-free embedding, is exactly zero.](figures/ctx_fig1.pdf){width=62%}

**Metric validity (anisotropy).** These representations are anisotropic, as expected [@ethayarajh2019]: the mean
cosine of random different-gene pairs is +0.45 to +0.53 at layers 4–11 (raw), and same-gene self-similarity
falls with depth (+0.94 → +0.82 → +0.75) — genes become more context-specific in upper layers, the single-cell
echo of Ethayarajh's finding. Our per-dimension z-scoring drives anisotropy to ~0.00 while keeping genes
distinguishable (corrected self-similarity +0.89/+0.65/+0.53), and reduces the top-PC share (0.21→0.12 at layer
11), so the cosine-based metrics are computed in a corrected space.

### 4.2 It is not the gene's rank position

MaxToki reads a cell as a rank-ordered gene list and encodes token position, so a gene ranking 5th in one cell
type and 500th in another is an obvious confound. It is not the driver (Fig. 2). The gene-specific response is
essentially identical for rank-stable genes (+0.759) and rank-moving genes (+0.755), and **96%** survives
regressing the projection on per-(gene, context) mean rank (+0.728 residualised at layer 4). Shift *magnitude*
correlates with rank change (ρ ≈ 0.35), but its *direction* is largely orthogonal to rank.

![**Figure 2.** The effect is not rank position: gene-specific response is identical for rank-stable and rank-moving genes and survives residualising on mean rank.](figures/ctx_fig2.pdf){width=62%}

### 4.3 It replicates on an independent dataset

All of the above is on one atlas. On an independent dataset — Setty CD34+ bone-marrow hematopoiesis (the same lineage in which a recoverable
manifold was found in scGPT [@kend2026manifold]), with eight *developmental states* rather than adult cell types as the context axis — the phenomenon replicates: EXCESS
**+0.640** (layer 2) and **+0.584** (layer 4), gene-shuffled null +0.004, functional-z **+9.1/+7.6**. Neither
the contextualisation nor its functional organisation is a Tabula-Sapiens or adult-cell-type artefact.

### 4.4 The modulation is functionally organised — but the organisation is co-expression

Context modulates genes along interpretable functional axes far more than along random directions. Reported as
raw z (rank-controlled in parentheses), the nuclear/transcriptional vs. surface/secreted axis reaches z =
**+21.9** (+20.2); mitochondrial vs. cytoskeletal +10.3 (+8.3); transcription vs. transport +7.1 (+6.3); the
abundance-neutral but weakly-separable cell-cycle-vs-differentiation axis is null, +0.7 (−0.7) at AUC 0.54. The
nuclear/surface axis is biologically valid — cell types order sensibly, lymphocytes at the nuclear pole (CD8 T
+4.1, CD4 T +3.7, B +3.4) and secretory/surface cells at the other (alveolar type-2 −5.3, type-1 −3.7, macrophage
−3.2).

The decisive control asks whether this exceeds **co-expression** (Fig. 3). Functional gene sets are co-regulated
modules, so their genes co-vary across cells regardless of any "understanding" of function. The right null must
match the functional poles on **both** their size and their co-expression: comparing a ~1,880-gene functional
axis against small random modules is invalid, because smaller sets reach higher coherence and different power.
We therefore build null gene-set pairs of the *same sizes* as the functional poles, matched on co-expression
coherence, and — because the null power distribution is heavily right-skewed — report the functional axis's power
as an empirical percentile among them. The nuclear/surface axis's power (5.32) sits at the **p = 0.05** border of
size-matched co-expression modules (mean 2.22); the mitochondrial/cytoskeletal axis at p = 0.16 and the
transcription/transport axis at p = 0.10 are not significant. So the functional organisation does *not* robustly
exceed co-expression — the headline axis is at best marginal, the others indistinguishable. It does clear a
representation-space **tightness**-matched null (mean cosine among the axis genes' own representations) and, by
construction, the random-axis null (p < 0.001), so it is not merely a compact-gene-set artefact — but
co-expression is sufficient. This is the geometry instance of the field-wide pattern that SCFMs re-express
co-expression [@kendiukhov2026attention; @kedzierska2025zeroshot].

![**Figure 3.** The crux control. Random gene modules (grey) trace context-modulation power against co-expression coherence; the functional axes (red) sit on that curve, so the functional organisation is co-expression.](figures/ctx_fig3.pdf){width=66%}

### 4.5 The functional-context direction is causally used

Level 1 is a genuine model-fact, and it is *used*, not decorative (Fig. 4). In a steer-here/read-there test we
add a functional direction *u* to a random half of a cell's genes and read the model's own output logits at the
*other*, untouched half. Nuclear-gene logits rise toward +*u* (specific effect **+0.160**) and fall toward −*u*
(**−0.112**), a signed swing of **+0.272** versus **−0.015** for a norm-matched random push, dose-dependently
across three strengths and beyond the random control on 24 of 30 cells (*p* = 7×10⁻⁴). The
functional-context direction propagates through attention and shifts the model's predictions.

The effect is not specific to one axis or layer. It replicates for the same nuclear/surface direction injected
at a deeper layer (layer 8: swing +0.566, 30/30 cells, *p* = 9×10⁻¹⁰, stronger than at layer 4) and for an
independent functional axis (mitochondrial vs. cytoskeletal at layer 4: swing +0.321, *p* = 0.02). It is null
for the weakest-separated axis we tested causally (transcription vs. transport, validity AUC 0.66: swing −0.226, *p* = 1.0 — in fact a reversal),
consistent with a channel that exists for well-formed functional directions but not for every ontology contrast
one can write down. Across the four causal configurations the p-values are raw; the two strongest (nuclear/surface at layers 4 and 8) survive Bonferroni correction for the four tests, while the mitochondrial axis (p = 0.02) becomes marginal.

![**Figure 4.** Causal use. Steering a functional direction into half a cell's genes shifts the model's logits at the *other* genes with the sign of the push, dose-dependently, far beyond a norm-matched random push.](figures/ctx_fig4.pdf){width=62%}

### 4.6 No novel context-dependent gene function (the co-expression ceiling)

Two independent attempts to read *context-appropriate or novel* gene function beyond co-expression are null.
(1) **Directional congruence** — does a gene move toward its own functional pole specifically in contexts sharing
that character (a cross-validated non-additivity test)? Inside the co-expression null (empirical p ≈ 0.2); no signal.
(2) **Placement near curated targets** — are transcription factors positioned near their ChIP-seq / literature
targets, context-appropriately (21 TFs)? Static excess +0.023 (*p* = 0.19), activity-modulation ρ = +0.080
(*p* = 0.33); no signal, even with real target ground truth. A third attempt (curated context-switching TFs)
is uninformative rather than null: only one such TF (RUNX1) is well-sampled in the panel, because context-
switching TFs are typically low-expressed — itself an instance of the measurement wall. The boundary is
consistent: an expression-only SCFM contextualises genes richly, but the information is co-expression, and does
not resolve a gene's changed role in a new cell beyond it (the ceiling is quantified in Fig. 3).

### 4.7 Contextualisation as a model metric: architecture, scale, and learned vs. architectural

The direction immune to the co-expression ceiling is to treat the phenomenon as a **model-quality metric**
(Table 2, Fig. 5). Two findings. (1) **Contextualisation is universal and scales.** Every architecture strongly
contextualises genes (EXCESS 0.70–0.88, main-effect replication > 0.93); on matched data (600 cells, identical
settings) it increases with size, MaxToki-217M +0.740 → 1B +0.881, in line with the loss-scaling behaviour reported for these masked-reconstruction models [@kend2026scaling]. (2) **Architecture determines functional
organisation.** Rank-based MaxToki organises context modulation along functional axes an order of magnitude more
(z = 16–22) than value-binned scGPT (+2.7) or protein-sequence-initialised STATE (+1.3). STATE is telling: it
moves genes strongly in context (EXCESS +0.835) but in no functionally coherent direction — its ESM2-initialised
gene tokens carry sequence identity, not co-expression. The gap is real, not low power (pole sizes 366–1,390).

(3) **The phenomenon is learned, not architectural.** An untrained, random-init MaxToki of identical
architecture and vocabulary — the decisive control, since attention mixes tokens regardless of training —
contextualises genes far less (EXCESS **+0.274** vs. +0.74–0.76 trained) and organises that modulation
functionally far less (functional-z **+3.8** vs. +21.2). This separates two things cleanly. Contextualisation has
a real architectural *floor* (+0.274: fixed random embeddings plus deterministic attention mixing give a
reproducible, gene-specific shift on their own), but roughly two-thirds of the trained effect is learned. The
functional organisation is *mostly* learned: the architecture alone produces a significant amount (z = +3.8),
which training then amplifies ~5.6× (to +21.2) — so roughly 18% of the functional organisation is architectural
and 82% is built by training, aligning context modulation with functional axes, and hence with co-expression.

| Model | tokenisation | contexts | genes | dim | EXCESS | functional-z |
|---|---|---|---|---|---|---|
| scGPT | value-binned | 8 | 5000 | 512 | +0.705 | +2.7 |
| STATE-SE (state-embedding) | ESM2 gene tokens | 9 | 1840 | 2048 | +0.835 | +1.3 |
| MaxToki-217M | rank-based | 12 | 5000 | 1232 | +0.740 | +21.2 |
| MaxToki-1B | rank-based | 12 | 5000 | 2304 | +0.881 | +16.4 |
| MaxToki-217M (random init) | untrained control | 12 | 6000 | 1232 | +0.274 | +3.8 |

*Table 2. Cross-model comparison at layer 4. EXCESS = contextualisation strength; functional-z = organisation
along the nuclear/surface axis vs. a random-axis null. MaxToki-217M/1B use matched 600-cell settings for a fair
scaling comparison; the random-init row is the learned-vs-architectural control (§4.7); its 6000-gene panel vs.
the 5000 of the trained rows does not affect EXCESS, which is a within-panel cosine statistic.*

![**Figure 5.** Contextualisation as a model metric across architecture, scale, and the random-weights control.](figures/ctx_fig5.pdf){width=92%}

### 4.8 Context aids prediction, but not through the most-contextualised genes

If contextualisation served the autoregressive objective directly, the genes the model contextualises most
should be the genes whose next-gene prediction most benefits from real context. It does not work that way. Using
a chimeric-context intervention — the model's log-probability of a target gene given its cell's true
rank-ordered predecessors versus a random other cell's — real context helps prediction enormously (mean benefit
**+7.7 nats**, a nat being a natural-log unit of log-probability; i.e. the true gene is roughly two-thousand-fold more probable under real than random context; 536
target genes). But the per-gene strength of that benefit does **not** track how much a gene is contextualised:
the partial Spearman correlation, controlling for gene frequency, is **−0.087** (95% CI [−0.168, −0.004]) — at
best faintly negative, and with per-gene estimates noisy on both sides, best read as *no positive link* rather
than a firm negative one. Context is central to the model's objective, and contextualisation is real, but the
two are not aligned per gene: contextualisation is not simply the model encoding its per-gene predictive reliance
on context.

## 5. Discussion

The result is a characterisation of *how* SCFMs represent genes in context, and a proof of *what* that
representation carries. Positively: SCFMs reshape gene representations gene-specifically, reproducibly, causally,
more with scale, and — for rank-based tokenisation — along co-expression-aligned functional axes. This is, to
our knowledge, the first measurement for single-cell models of the "how contextual are contextualised
representations" question [@ethayarajh2019], and the architecture/scale comparison is a usable interpretability
metric that discriminates tokenisation schemes a downstream benchmark would not.

Negatively, and equally a contribution: the modulation does not carry a gene's context-specific function beyond
co-expression. The exciting version — reading a gene's changed role in a new cell to generate hypotheses — is
not reachable for an expression-only model by either route we tried, each with a matched-null control. Why the
ceiling exists is mechanistically coherent: the model's only inputs are which genes are co-expressed, so the
richest thing its contextual geometry can encode about a gene is its co-expression neighbourhood; genuinely
novel context-dependent function is frequently post-transcriptional [@Jeffery1999Moonlighting; @Heinz2010LDTFs]
and leaves no expression signature to latch onto. The finding therefore predicts where headroom is: multi-omic
or sequence-aware models, not larger expression-only ones, are what a Level-2 result would require.

## 6. Limitations

Cross-architecture EXCESS is confounded by occurrence cap (20 for scGPT/STATE vs. 50 for MaxToki) and
dimensionality (higher-dimensional models show higher EXCESS), so the 217M→1B increment reflects size *and*
width, and even the cap-matched scaling pair still differs in width (dim 1232 vs. 2304) — the width-robust
comparison is the functional-z architecture gap, which has ample power. The functional-z gap itself is measured
on non-identical context sets across models (12 vs. 9 vs. 8 cell types); we did not re-run it on a common context
set, so a residual context-composition confound remains. The causal test (§4.5) compares the functional
direction only against a norm-matched *random* push, not against co-expression-coherence-matched directions, so
it shows the direction is used, not that its causal use exceeds co-expression. STATE is evaluated on immune
cells only in this study, so its cross-architecture point rests on a narrower panel than the others. The co-expression ceiling is a null result: the headline axis sits
at the p = 0.05 border of the size-matched co-expression null and the other axes are non-significant, which shows
the organisation does not clearly exceed co-expression, not that the two are provably identical; the empirical
p depends on the null-module construction. The
Level-2 null is evidence of a ceiling, not proof of absence — some context-dependent function is structurally
invisible to an expression-only model, so those tests are weak by construction (pre-registered). Functional axes
are built from one ontology file; the verdict is robust across three axes but specific z-scores depend on pole
definitions.

## 7. Conclusions

Single-cell foundation models build contextual gene representations that are a genuine object of study: a gene's
vector is reshaped by its cellular context in a gene-specific, reproducible, learned, and causally used way that
grows with model scale and is shaped by tokenisation. This is the first systematic characterisation of the
phenomenon for single-cell models, and it yields a usable model metric — contextualisation strength and its
functional organisation — that distinguishes architectures a downstream benchmark would not. At the same time,
matched-null controls place a firm boundary on what that contextualisation carries: it does not robustly exceed
co-expression, and two independent probes for novel context-dependent gene function are null. The practical
implication is that reading a gene's changed role in a new cell type — the application that would make these
representations a hypothesis-generation tool — is not reachable with an expression-only model, and will require
multi-omic or sequence-aware models rather than merely larger expression-trained ones.

## 8. Methods details and reproducibility

All analyses are deterministic (seed 0) and run on cached extractions. Extraction: `ctx_extract_maxtoki.py`
(+ `_scgpt`, `_state`, `_devel`, `_random`). Analyses: `ctx_polysemy.py` (§4.1), `ctx_anisotropy.py` (§4.1),
`ctx_layer_curve.py` (§4.1), `ctx_position_confound.py` (§4.2), `ctx_independent.py` (§4.3),
`ctx_functional_axes.py`, `ctx_coexpr_null_v2.py`, `ctx_tightness_null.py` (§4.4), `ctx_causal.py` across
axes/layers (§4.5), `ctx_directional_probe.py` / `ctx_curated_targets.py` (§4.6), `ctx_cross_model.py` (§4.7),
`ctx_prediction_link.py` (§4.8). Metric definitions: EXCESS as in §3; reproducible interaction power = mean over
count-balanced (gene, context) of the product of the two partitions' two-way interaction terms along a unit
axis; functional-z, its z-score against 200–300 random-partition axes; co-expression coherence = mean within-set
expression correlation. The co-expression null (`ctx_coexpr_null_v2.py`) draws null gene-set pairs of the *same
sizes* as the functional poles, matched on coherence, and reports the functional axis's power as an empirical
percentile among them (skew-robust). Result artefacts are the corresponding `results/ctx_*.json`.

## CRediT authorship contribution statement

**Ihor Kendiukhov:** Conceptualization, Methodology, Software, Formal analysis, Investigation, Data curation,
Writing – original draft, Writing – review & editing, Visualization.

## Declaration of competing interest

The author declares that he has no known competing financial interests or personal relationships that could have
appeared to influence the work reported in this paper.

## Data and code availability

All analysis code is openly available at <https://github.com/Biodyn-AI/gene-context>, including every analysis
script, the intermediate result artefacts (`results/ctx_*.json`), the figure-generation code, and this
manuscript with its bibliography. The study uses only publicly available models and data: the MaxToki-217M/1B
checkpoints (Hugging Face, `theodoris-lab/MaxToki`); scGPT, Geneformer and STATE (public releases); the Tabula
Sapiens atlas [@tabulasapiens2022] (CZ CELLxGENE); the Setty CD34+ bone-marrow dataset; and Gene Ontology
annotations. The per-gene contextual representation tensors extracted from these models (the intermediate
`.npz` files, several gigabytes) are not deposited because of their size, but are fully regenerable from the
public inputs with the provided extraction scripts, and are available from the author on reasonable request. No
new experimental data were generated.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or
not-for-profit sectors.

## References
