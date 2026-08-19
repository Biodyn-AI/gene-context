# Cover letter

Dear Editor,

I am pleased to submit the manuscript **"Single-cell foundation models represent genes in context, but the
information is co-expression"** for consideration as a research article in *Computational Biology and Chemistry*.

Single-cell foundation models (SCFMs) such as Geneformer, scGPT and the recent MaxToki are increasingly used to
produce gene and cell representations, yet what their *contextual* gene representations — a gene's vector inside
a particular cell — actually encode has not been characterised. This manuscript provides that characterisation,
and does so with the matched-null rigour that the field's recent evaluation literature has shown to be
necessary.

The main contributions are:

1. A systematic demonstration that an SCFM reshapes a gene's representation by its cellular context in a way that
   is gene-specific, reproducible across independent cells, not an artefact of token rank, replicable on an
   independent dataset, causally used by the model, learned rather than architectural, and increasing with model
   scale — the single-cell analogue of the "how contextual are contextualised representations" question from
   language-model interpretability.

2. A usable model-characterisation metric: contextualisation strength and its functional organisation differ by
   an order of magnitude across tokenisation schemes (rank-based vs. value-binned vs. protein-sequence-token),
   distinguishing architectures in a way a downstream task benchmark does not.

3. A carefully controlled negative result: against size-matched co-expression nulls, the functional organisation
   of that contextualisation does not robustly exceed co-expression, and two independent probes for *novel*
   context-dependent gene function are null. This places a precise, reproducible boundary on what expression-only
   SCFMs can offer as hypothesis-generation tools, and predicts where the headroom lies (multi-omic or
   sequence-aware models).

We believe the work fits the scope of *Computational Biology and Chemistry* at the intersection of machine
learning and computational biology, and will be of interest to readers building or evaluating single-cell
foundation models. The study uses only publicly available models and data; all code and intermediate result
artefacts are openly available at <https://github.com/Biodyn-AI/gene-context>.

This manuscript is original, has not been published previously, and is not under consideration for publication
elsewhere. The author declares no competing interests.

Thank you for your consideration.

Sincerely,

Ihor Kendiukhov
