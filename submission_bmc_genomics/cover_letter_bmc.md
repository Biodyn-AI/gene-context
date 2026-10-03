Ihor Kendiukhov\
Institute of Medical Genetics and Applied Genomics, University of Tübingen\
Calwerstraße 7, 72076 Tübingen, Germany\
kendiukhov@gmail.com\
ORCID: 0000-0001-5342-1499

[Date]

The Editor\
*BMC Genomics*

Dear Editor,

I am pleased to submit the manuscript **"Single-cell foundation models represent genes in context, but reveal no context-specific gene function beyond co-expression"** for consideration as a Research article in *BMC Genomics*. I suggest the **Transcriptomic Methods** section.

Single-cell foundation models such as Geneformer, scGPT and the recent MaxToki are now widely used to produce gene and cell representations from single-cell RNA-sequencing data. Inside these models, every gene has a different vector in every cell. Whether that contextual representation shows how a gene's role differs between cell types, beyond what the expression data already show, has not been tested. This manuscript tests it, with the matched controls that recent evaluations of these models have shown to be necessary.

The main results are:

1. **The contextual representation is real, in every model tested, but not unique to models.** MaxToki-217M, MaxToki-1B, scGPT and STATE, run on identical cells, all change a gene's representation with cell type in a gene-specific way that reproduces across independent cells (agreement between independent halves of the cells 0.74 to 0.87 above a shuffled control) and is not explained by the gene's rank in the cell. In MaxToki-217M it replicates on an independent bone-marrow dataset. Untrained models keep 36–88% of it, and a representation built from expression alone, with no model, also shows it (+0.321 when each gene vector averages 50 cells, as in the models).
2. **What the models add is organisation by co-expression.** In MaxToki the change follows functional axes (for example nuclear versus surface genes); STATE shows weaker organisation on that axis and scGPT little. In no model is this organisation significantly stronger than that of gene sets matched on size and co-expression. In every model the change follows co-expressed gene sets; in three representations built from expression alone, on the same cells, it does not. Steering along the nuclear/surface direction changes MaxToki's predictions for untouched genes in both sizes, more than directions built from random groups of the same genes; MaxToki-1B also reads the other two functional directions tested, while MaxToki-217M reads one of them in reverse. The pushes are larger than natural changes, so this shows that the models can read these directions, not that they use the natural changes. In both MaxToki sizes, genes that change more gain slightly more from their own context in next-gene prediction.
3. **No context-specific gene function beyond co-expression.** Two direct tests, run in all four models, are negative, although they could detect only clear effects. The only signal beyond co-expression is a static closeness of transcription factors to their curated targets that does not change with cell type; in STATE it is already present in the untrained model, which shares the trained model's protein-sequence gene table. The functional organisation of contextual change separates some architectures but does not yet rank them reliably.

I believe the work fits *BMC Genomics* and its Transcriptomic Methods section: it is a methods-focused, fully reproducible evaluation of tools that are increasingly used to analyse transcriptomic data, and it reports both positive and negative results with explicit controls. All code, the result files from which every number, table and figure is computed, and the figure code are openly available at <https://github.com/Biodyn-AI/gene-context>.

**Editorial policies.** The study is a secondary analysis of publicly available, de-identified human single-cell data and public models; no new human data were collected, so ethics approval was not required. The use of an AI assistant in coding, reference checking and drafting is described in the Methods, as BMC policy requires. Related work by the author on static gene embeddings, attention, residual-stream geometry, the geometry and topology of hidden states and sparse-autoencoder features is cited in the manuscript (six papers). Three of them are published: in *BMC Genomics* (doi:10.1186/s12864-026-12965-8), *BMC Bioinformatics* (doi:10.1186/s12859-026-06538-5) and *PLOS One* (doi:10.1371/journal.pone.0344826); the other three are preprints, and the sparse-autoencoder atlas paper also has a Research Square preprint (doi:10.21203/rs.3.rs-9082479/v1). None of them reports the contextual-representation analyses presented here. [Author: state the current status of the sparse-autoencoder atlas paper, and confirm that this manuscript has not been posted as a preprint and that the earlier Elsevier submission was not made or has been withdrawn.]

**Competing interests.** I declare no competing interests.

I confirm that I am the sole author, that I have read and approved the manuscript for submission, and that the manuscript has not been published and is not under consideration for publication elsewhere.

Thank you for your consideration.

Sincerely,

Ihor Kendiukhov
