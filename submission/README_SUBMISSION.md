# Submission package — Computational Biology and Chemistry (Elsevier)

Manuscript: *Single-cell foundation models represent genes in context, but reveal no context-specific gene function
beyond co-expression* (I. Kendiukhov). Prepared under Elsevier's "Your Paper Your Way" (free-format initial
submission). Rebuilt on 3 Oct 2026 from the same text and numbers as the BMC Genomics package.

**This is an alternative to `../bmc_submission/`, not an addition.** Both packages describe the same study. Submit
only one of them at a time; the cover letters say the manuscript is not under consideration elsewhere.

## How it is built

`bash cbc_submission/build/build_cbc.sh` (run from `route_genemanifold/`) takes the rendered BMC manuscript
(`../bmc_submission/manuscript_bmc.md`, itself built from the result files) and changes only what Elsevier needs:
a one-paragraph abstract of at most 250 words (`cbc_abstract.template.md`), 3–5 highlights of at most 85 characters
(`highlights.template.txt`), 7 keywords, "Introduction" instead of "Background", Elsevier end statements (CRediT,
competing interests, generative-AI declaration, data availability, funding) in place of BMC's Declarations, the
supplement as "Supplementary Table S1", and figures placed with their captions for the review PDF. The build fails if
a limit is broken. So rebuild the BMC package first if any result changes.

## Files and what to upload where

| File | Editorial Manager item |
|---|---|
| `manuscript.pdf` | Manuscript (review PDF: continuous line numbers, figures with captions, references) |
| `manuscript.tex` | LaTeX source (optional at initial submission; required at revision) |
| `title_page.docx` / `title_page.pdf` | Title page |
| `cover_letter.docx` / `cover_letter.pdf` | Cover letter |
| `highlights.txt` | Highlights (5 bullets, each ≤ 85 characters, checked at build) |
| `title_abstract_keywords.txt` | Title / abstract / keywords for the submission form (one per line) |
| `declaration_of_interest.txt` | Declaration of interests |
| `credit_author_statement.txt` | CRediT authorship statement |
| `figures/Figure1..5.pdf` | Separate figure files (vector PDF) |
| `Supplementary_Table_S1.xlsx` | Supplementary material (the full multi-model results; same file as the BMC Additional file 1) |

## Before submitting — fill in or confirm

1. The date in the cover letter, and your ORCID in the submission system.
2. The bracketed instruction at the end of the cover letter's editorial-policies paragraph (status of the
   sparse-autoencoder atlas paper; confirm this manuscript is not a preprint and is not under consideration at BMC
   Genomics). Then delete the bracket.
3. The archived-software DOI: the data availability statement still says "DOI [to be added before submission]".
   Make a Zenodo release of the tag `bmc-submission-v3.1` of `Biodyn-AI/gene-context` and put the DOI in the BMC
   template, then rebuild both packages.
4. The generative-AI declaration and the funding / competing-interest wording.

See `../bmc_submission/README_SUBMISSION.md` for what changed since the 1 Oct version.
