# Submission package — BMC Genomics

Manuscript: *Single-cell foundation models represent genes in context, but reveal no context-specific gene
function beyond co-expression* (I. Kendiukhov). Article type: **Research article**. Suggested section:
**Transcriptomic Methods**.

Prepared against the BMC Genomics submission guidelines (research article) as of 1 Oct 2026.

## Files and where they go in the submission system

| File | Upload as |
|---|---|
| `manuscript_bmc.docx` | **Manuscript** (main file; editable DOCX as BMC requires). Title page, abstract, all sections, declarations, references, table titles/legends and figure titles/legends are inside it. Double line spacing, continuous line numbers, page numbers, no page breaks. |
| `figures/Figure1.pdf` … `figures/Figure5.pdf` | **Figures**, one file each (vector PDF, TrueType fonts embedded, no titles inside the graphic; 85 mm or 170 mm wide). |
| `cover_letter_bmc.docx` / `cover_letter_bmc.pdf` | **Cover letter** (paste or upload; same text). |
| `title_abstract_keywords.txt` | Text for the online form: line 1 title, line 2 abstract, line 3 keywords. |
| `manuscript_bmc_preview.pdf` | For your own checking only — a PDF of the same text. Do **not** upload it as the manuscript (BMC needs the editable file). |
| `manuscript_bmc.template.md`, `cover_letter_bmc.template.md`, `references_bmc.bib`, `build/` | Sources. Every number in the text is a `{{key}}` placeholder filled from `results/ctx_*.json` by `build/make_numbers.py`. `bash build/build_bmc.sh` regenerates the numbers, renders `manuscript_bmc.md` and `cover_letter_bmc.md`, redraws the figures and builds the DOCX, preview PDF and cover letter. |

There are no additional files and no tables outside the manuscript (Tables 1–3 are each under one page and sit
in the text).

## BMC requirements checked

- [x] Main file editable (DOCX). Not LaTeX, because BMC compiles LaTeX with pdfLaTeX and our builds use XeLaTeX.
- [x] Title page: title, author full name, institution, corresponding author email.
- [x] Abstract ≤ 350 words with Background / Results / Conclusions; no references (checked at every build).
- [x] 3–10 keywords (8).
- [x] Section order: Background, Methods, Results, Discussion (with Limitations), Conclusions, List of
      abbreviations, Declarations, References, Figure titles and legends.
- [x] All seven required declaration headings present (optional "Authors' information" omitted); "Not applicable" where relevant.
- [x] Ethics statement for public, de-identified human data.
- [x] Use of a large language model documented in Methods ("Use of large language models").
- [x] Availability of data and materials in BMC form, with every dataset and model cited in the reference list
      with a full persistent link; software availability block (project name, home page, archived version,
      OS, language, requirements, licence, restrictions).
- [x] References: numbered Vancouver style in square brackets (BMC CSL); every web link is a numbered reference
      with an access date; every DOI checked against CrossRef/DataCite on 1 Oct 2026.
- [x] Tables: title above (≤ 15 words), legend below, no shading, plain black rules, no commas in numbers.
- [x] Figures: separate files, fonts embedded, width ≤ 170 mm, titles ≤ 15 words and legends ≤ 300 words in the
      manuscript.
- [x] Double spacing, line numbers, page numbers, no page breaks.

## Before you submit — you must fill in or confirm

1. **The date** in the cover letter. Address and email are now the Tübingen affiliation and kendiukhov@gmail.com (as
   in your Scientific Reports package); add your ORCID in the submission system.
2. **Preprint / related submissions**: the cover letter ends its editorial-policies paragraph with one bracketed
   instruction: state the status of your sparse-autoencoder atlas paper (Research Square preprint, "In Review"), and
   confirm this manuscript is not a preprint and the Elsevier submission was not made or was withdrawn. Then delete the
   bracket.
3. **Not submitted elsewhere.** The cover letter states the manuscript is not under consideration elsewhere. A separate
   package for *Computational Biology and Chemistry* (Elsevier) was prepared earlier and is now superseded. Submit to BMC
   only if that submission was not made or has been withdrawn.
4. **LLM statement** (Methods, "Use of large language models") — confirm the wording matches how the assistant was used.
5. **Archived software version**: the manuscript names the git tag `bmc-submission-v1`. For a DOI, enable the
   Zenodo–GitHub integration for `Biodyn-AI/gene-context`, make a release from that tag, and replace the line with the DOI.
6. **Competing interests / funding**: confirm "none" and "Not applicable".
7. Optional: suggested reviewers (the submission system asks; BMC does not require them).

## What changed relative to the earlier (Elsevier) version of the paper

Found while preparing this package, by checking the code and by two independent agent reviews. All are fixed here and in
the general manuscript in the repository.

- **Input encoding (the big one).** The Tabula Sapiens files store log1p(counts per 10,000) in `X`; the old code treated
  that as counts and log-transformed it a second time before ranking, so MaxToki saw a gene order driven mostly by global
  gene medians. All MaxToki inputs now use the model's own encoding (raw counts / gene median, `ctx_tokenise.py`,
  checked against the MaxToki tokenizer: identical up to the order of exactly tied genes), and every MaxToki extraction
  and analysis was re-run (`rerun_tokfix.sh`, 1 Oct 2026; the spliced-only Setty robustness steps at its end were run
  separately on 2 Oct). Several results changed:
  - functional organisation is now inside both co-expression nulls (matched p = 0.08–0.38; strong p = 0.70–0.80);
  - the prediction link reversed: genes that change more with context gain more from real context (partial ρ +0.16);
  - steering works for the nuclear/surface axis (mitochondrion/cytoskeleton: no consistent effect, significant only at the
    largest strength, p = 0.049 uncorrected; transcription/transport: reversed);
  - no consistent scaling from 217M to 1B once relative depth is matched (the old "0.74 → 0.88" was layer 4 vs layer 4);
  - TFs sit slightly closer to curated targets than co-expression predicts (static), but not context-appropriately.
- **Other method fixes:** null gene sets matched on size and co-expression (and on tightness) with disjoint poles and
  Monte Carlo p values; gene-level bootstrap CIs; depth-matched scaling; random-weights control on the same cells and
  panel as its comparator; the full 5,780-cell Setty data with spliced + unspliced counts (spliced-only as robustness);
  a cleaner chimeric control for the prediction link; deterministic rank control.
- **Text fixes:** layer peak (layers 1–2), STATE's cell types (immune, kidney, lung), scGPT and depth caveats, which
  analyses were run on which model, reference updates (journal versions of your attention paper, UCE and STATE; your
  BMC Bioinformatics paper added; versioned dataset links).
