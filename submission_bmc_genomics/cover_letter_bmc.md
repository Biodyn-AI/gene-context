Ihor Kendiukhov\
Institute of Medical Genetics and Applied Genomics, University of Tübingen\
Calwerstraße 7, 72076 Tübingen, Germany\
kendiukhov@gmail.com

[Date]

The Editor\
*BMC Genomics*

Dear Editor,

I am pleased to submit the manuscript **"Single-cell foundation models represent genes in context, but reveal no context-specific gene function beyond co-expression"** for consideration as a Research article in *BMC Genomics*. I suggest the **Transcriptomic Methods** section.

Single-cell foundation models such as Geneformer, scGPT and the recent MaxToki are now widely used to produce gene and cell representations from single-cell RNA-sequencing data. Inside these models, every gene has a different vector in every cell. Whether that contextual representation tells us anything about a gene's role in a particular cell type, beyond what the expression data already show, has not been tested. This manuscript tests it, with the matched controls that recent evaluations of these models have shown to be necessary.

The main results are:

1. **The contextual representation is real.** In MaxToki-217M, a gene's representation changes with cell type in a gene-specific way that reproduces across independent cells (cosine 0.76 against 0.001 for a shuffled control), is strongest in the first two layers, is not explained by the gene's rank in the cell, and replicates on an independent bone-marrow dataset. Training strengthens it (an untrained model of the same architecture keeps about 60% of the effect but only about a third of its functional organisation), and it is weakly aligned with the model's objective: genes that change more with context gain more from real context in next-gene prediction. All four models tested (MaxToki-217M and -1B, scGPT, STATE) contextualise genes strongly.
2. **Nothing context-specific beyond co-expression.** The context changes follow functional axes, but gene sets matched on size and co-expression give similar organisation (p = 0.08–0.38). Two direct tests for context-specific gene function, run in MaxToki-217M, are negative; the only signal beyond co-expression is a small, static closeness of transcription factors to their curated targets. Steering shows that the model reads one functional direction (nuclear versus surface genes); a second direction has no consistent effect and a third acts in reverse. Reading a gene's changed role in a new cell type is therefore not reachable with the model tested.
3. **A model metric.** The functional organisation of contextual change differs sharply between architectures (z = 8.1–18.4 for MaxToki, 3.2 for scGPT, 1.1 for STATE), while all models contextualise genes strongly.

I believe the work fits *BMC Genomics* and its Transcriptomic Methods section: it is a methods-focused, fully reproducible evaluation of tools that are increasingly used to analyse transcriptomic data, and it reports both positive and negative results with explicit controls. All code, the result files from which every number, table and figure is computed, and the figure code are openly available at <https://github.com/Biodyn-AI/gene-context>.

**Editorial policies.** The study is a secondary analysis of publicly available, de-identified human single-cell data and public models; no new human data were collected, so ethics approval was not required. The use of an AI assistant in coding, reference checking and drafting is described in the Methods, as BMC policy requires. Related work by the author on static gene embeddings, attention, residual-stream geometry, the geometry and topology of hidden states, a haematopoietic manifold in scGPT, sparse-autoencoder features and scaling laws is cited in the manuscript (eight papers). Four of them are published: in *BMC Genomics* (doi:10.1186/s12864-026-12965-8), *BMC Bioinformatics* (doi:10.1186/s12859-026-06538-5), *PLOS One* (doi:10.1371/journal.pone.0344826) and *Transactions on Machine Learning Research*; the other four are preprints, and the sparse-autoencoder atlas paper also has a Research Square preprint (doi:10.21203/rs.3.rs-9082479/v1). None of them reports the contextual-representation analyses presented here. [Author: state the current status of the sparse-autoencoder atlas paper, and confirm that this manuscript has not been posted as a preprint and that the earlier Elsevier submission was not made or has been withdrawn.]

**Competing interests.** I declare no competing interests.

I confirm that I am the sole author, that I have read and approved the manuscript for submission, and that the manuscript has not been published and is not under consideration for publication elsewhere.

Thank you for your consideration.

Sincerely,

Ihor Kendiukhov
