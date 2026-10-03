# Submission package — BMC Genomics

Manuscript: *Single-cell foundation models represent genes in context, but reveal no context-specific gene
function beyond co-expression* (I. Kendiukhov). Article type: **Research article**. Suggested section:
**Transcriptomic Methods**.

Prepared against the BMC Genomics submission guidelines (research article) on 1 Oct 2026; rebuilt on 3 Oct 2026
after the multi-model analyses and three rounds of claim checking (see "What changed" below).

## Files and where they go in the submission system

| File | Upload as |
|---|---|
| `manuscript_bmc.docx` | **Manuscript** (main file; editable DOCX as BMC requires). Title page, abstract, all sections, declarations, references, Tables 1–5 with their titles and legends, and figure titles/legends are inside it. Double line spacing, continuous line numbers, page numbers, no page breaks. |
| `figures/Figure1.pdf` … `figures/Figure5.pdf` | **Figures**, one file each (vector PDF, fonts embedded, no titles inside the graphic; 85 mm or 170 mm wide). |
| `Additional_file_1.xlsx` | **Additional file 1** (Excel workbook): the full results of the multi-model comparison, one sheet per analysis, generated from the result files by `build/make_additional_file.py`. Its title and description are in the manuscript under "Additional files". |
| `cover_letter_bmc.docx` / `cover_letter_bmc.pdf` | **Cover letter** (paste or upload; same text). |
| `title_abstract_keywords.txt` | Text for the online form: line 1 title, line 2 abstract, line 3 keywords. |
| `manuscript_bmc_preview.pdf` | For your own checking only — a PDF of the same text. Do **not** upload it as the manuscript (BMC needs the editable file). |
| `manuscript_bmc.template.md`, `cover_letter_bmc.template.md`, `references_bmc.bib`, `build/` | Sources. Every number in the text is a `{{key}}` placeholder filled from `results/ctx_*.json` by `build/make_numbers.py`. `bash build/build_bmc.sh` regenerates the numbers, renders `manuscript_bmc.md` and `cover_letter_bmc.md`, redraws the figures, writes Additional file 1 and builds the DOCX, preview PDF and cover letter. |

## BMC requirements checked

- [x] Main file editable (DOCX). Not LaTeX, because BMC compiles LaTeX with pdfLaTeX and our builds use XeLaTeX.
- [x] Title page: title, author full name, institution, corresponding author email.
- [x] Abstract ≤ 350 words with Background / Results / Conclusions; no references (checked at every build).
- [x] 3–10 keywords (8).
- [x] Section order: Background, Methods, Results, Discussion (with Limitations), Conclusions, List of
      abbreviations, Declarations, References, Figure titles and legends, Additional files.
- [x] All seven required declaration headings present (optional "Authors' information" omitted); "Not applicable" where relevant.
- [x] Ethics statement for public, de-identified human data.
- [x] Use of a large language model documented in Methods ("Use of large language models").
- [x] Availability of data and materials in BMC form, with every dataset and model cited in the reference list
      with a full persistent link; software availability block (project name, home page, archived version,
      OS, language, requirements, licence, restrictions).
- [x] References: numbered Vancouver style in square brackets (BMC CSL); every web link is a numbered reference
      with an access date.
- [x] Tables cited in order (1–5), title above (≤ 15 words), legend below, no shading, plain rules.
- [x] Figures: separate files, fonts embedded, width ≤ 170 mm, titles ≤ 15 words and legends ≤ 300 words in the
      manuscript.
- [x] Additional file: named "Additional file 1", with title and description in the manuscript, cited in the text.
- [x] Double spacing, line numbers, page numbers, no page breaks.

## Before you submit — you must fill in or confirm

1. **The date** in the cover letter. Your ORCID (0000-0001-5342-1499) is now in the cover letter; also enter it in the submission system.
2. **Preprint / related submissions: done.** The cover letter names the six related papers cited in the manuscript
   (three published, three preprints; the sparse-autoencoder atlas paper's Research Square entry is a plain posted
   preprint, not under review at a journal), links your Google Scholar profile for the rest, and states that this
   manuscript has not been posted as a preprint. If you opt into BMC's "In Review" service during submission, the
   manuscript will be posted on Research Square; that is fine, because the letter describes its status at submission.
3. **One journal only.** You did not submit to *Computational Biology and Chemistry*, so nothing needs withdrawing.
   The package in `../cbc_submission/` stays as an alternative; do not submit both.
4. **LLM statement** (Methods, "Use of large language models") — confirm the wording matches how the assistant was used.
5. **Archived software version: done.** Zenodo archived the GitHub release `release2` (commit `509f4b9`, the same
   code as tag `bmc-submission-v3.1`). The manuscript cites the version DOI 10.5281/zenodo.23119961 in the
   Software block, the data availability statement and the reference list; 10.5281/zenodo.23119960 is the
   all-versions DOI used in the repository badge. Do not make further GitHub releases before submission: each one
   creates a new Zenodo version.
6. **Competing interests / funding**: confirm "none" and "Not applicable".
7. Optional: suggested reviewers (the submission system asks; BMC does not require them).

## What changed since the 1 Oct package

- **Four models on the same cells.** scGPT and STATE were re-extracted with their standard inputs (scGPT: binned
  counts of all expressed genes with `<cls>`; STATE: its own data collator), on exactly the same 600 cells per cell
  type, partitions, gene panel and cap as MaxToki-217M and MaxToki-1B. Every Level 1 and Level 2 test (rank
  control, all four nulls, directional congruence, curated targets) now runs on all four models (Tables 3 and 4),
  with untrained versions of three architectures (Table 5).
- **Direct expression-only test.** Three representations built from expression alone, with no model, go through the
  identical pipeline. The centroid version was first built with 512 principal components; it failed the positive
  control (the shared shift did not replicate between cell halves), so it was rebuilt with 50 components and with at
  most 50 cells per gene, the same sample size as the models' 50 occurrences. Built this way, expression alone gives a
  gene-specific change (EXCESS +0.32, against +0.74 for MaxToki-217M) but no functional organisation, and co-expressed
  gene sets move no more than random ones.
- **Steering** now has a second control matched to the axis (30 fixed random splits of the pole genes), a push unit
  that is fair across sizes, MaxToki-1B at matched depth, and a check with the two very large hidden dimensions
  removed.
- **Prediction link** now uses a paired design (same cells, positions and donors for both MaxToki sizes) and a
  same-cell-type donor arm.
- **Three rounds of claim checking** against the result files (100, 74 and 23 confirmed issues) led to many wording changes:
  corrected multiple-testing families, the TF-specificity results reported with their corrections, caveats on the
  steering push sizes and on what the tests do not show, tissue and donor confounds stated with numbers, and the
  depth pattern of each model.
- **Additional file 1** (new) holds the full tables behind every multi-model number.

## Earlier changes (relative to the Elsevier version of 30 Sep)

- **Input encoding.** The Tabula Sapiens files store log1p(counts per 10,000) in `X`; the old code treated that as
  counts and log-transformed it a second time before ranking. All MaxToki inputs now use the model's own encoding
  (raw counts / gene median, `ctx_tokenise.py`), and every MaxToki extraction and analysis was re-run.
- Null gene sets matched on size and co-expression (and on tightness) with disjoint poles and Monte Carlo p values;
  gene-level bootstrap CIs; depth-matched scaling; random-weights control on the same cells and panel; the full
  Setty data with spliced + unspliced counts; deterministic rank control.
