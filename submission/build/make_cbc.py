"""Build the Computational Biology and Chemistry (Elsevier, "Your Paper Your Way") manuscript from the rendered BMC
manuscript, so both versions share one source of text and numbers.

Differences from the BMC version: one-paragraph abstract (<= 250 words, from cbc_abstract.template.md), up to 7
keywords, "Introduction" instead of "Background", Elsevier end statements (CRediT, competing interest, generative-AI
declaration, data availability, funding) in place of BMC's Declarations, no abbreviation list, and figures placed with
their captions for the review PDF. Writes cbc_submission/manuscript_cbc.md, highlights.txt, title_abstract_keywords.txt,
credit_author_statement.txt, declaration_of_interest.txt and title_page.md. Fails if a limit is broken."""
import os, re, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); CBC = os.path.dirname(HERE); RG = os.path.dirname(CBC)
BMC = os.path.join(RG, "bmc_submission")
if not os.path.isdir(BMC):                       # public repository layout
    BMC = os.path.join(RG, "submission_bmc_genomics")
N = json.load(open(os.path.join(BMC, "build", "numbers.json")))


def fill(path):
    s = open(path).read()
    miss = sorted({k for k in re.findall(r"\{\{([A-Za-z0-9_]+)\}\}", s) if k not in N})
    if miss:
        sys.exit(f"unresolved placeholders in {path}: {miss}")
    return re.sub(r"\{\{([A-Za-z0-9_]+)\}\}", lambda m: str(N[m.group(1)]), s)


s = open(os.path.join(BMC, "manuscript_bmc.md")).read()
title = re.search(r'title: "(.*)"', s).group(1)
aff = re.search(r"\^1\^ (.*)\n", s).group(1)
corr = re.search(r"Correspondence: (\S+)", s).group(1)

# ---- abstract, keywords, highlights ----
abstract = fill(os.path.join(CBC, "cbc_abstract.template.md")).strip()
n_words = len(re.sub(r"[*^]", "", abstract).split())
if n_words > 250:
    sys.exit(f"CBC abstract has {n_words} words (limit 250)")
kw_all = [k.strip() for k in s[s.index("**Keywords:**") + 13:].split("\n")[0].split(";")]
KEEP = ["Single-cell foundation models", "Contextual gene representations", "Mechanistic interpretability",
        "Gene co-expression", "Gene embeddings", "Transformers", "Model evaluation"]
kw = [k for k in KEEP if k in kw_all]
assert 1 <= len(kw) <= 7, kw
hl = [l.strip() for l in fill(os.path.join(CBC, "highlights.template.txt")).splitlines() if l.strip()]
bad = [h for h in hl if len(h) > 85]
if not 3 <= len(hl) <= 5 or bad:
    sys.exit(f"highlights must be 3-5 lines of <= 85 characters; too long: {bad}")

# ---- body: from Background to before Declarations; figures with captions ----
body = s[s.index("# Background"):s.index("# List of abbreviations")]
body = body.replace("# Background", "# Introduction", 1)
# Elsevier names supplements "Supplementary material"; the BMC Additional file 1 ships as Supplementary Table S1
body = body.replace("Additional file 1", "Supplementary Table S1")
supp = s[s.index("**Additional file 1**"):].strip().replace("**Additional file 1** (Excel workbook, .xlsx).",
                                                                 "**Supplementary Table S1** (Excel workbook, Supplementary_Table_S1.xlsx).")
import shutil
shutil.copy(os.path.join(BMC, "Additional_file_1.xlsx"), os.path.join(CBC, "Supplementary_Table_S1.xlsx"))
decl = s[s.index("# Declarations"):s.index("# References")]
avail = decl[decl.index("## Availability of data and materials") + len("## Availability of data and materials"):
             decl.index("## Competing interests")].strip()
ack = decl[decl.index("## Acknowledgements") + len("## Acknowledgements"):].strip()
ethics = decl[decl.index("## Ethics approval and consent to participate") + 45:decl.index("## Consent for publication")].strip()
avail = avail.replace("git tag `bmc-submission-v1`", "the git tag of the submitted version")
end = f"""# CRediT authorship contribution statement

**Ihor Kendiukhov:** Conceptualization, Methodology, Software, Formal analysis, Investigation, Data curation, Visualization, Writing – original draft, Writing – review & editing.

# Declaration of competing interest

The author declares that he has no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

# Declaration of generative AI and AI-assisted technologies in the writing process

During the preparation of this work the author used Claude (Anthropic) to help write and run the analysis code, check references, and draft and edit the text. After using this tool, the author reviewed and edited the content as needed and takes full responsibility for the content of the published article.

# Ethics statement

{ethics}

# Data availability

{avail}

# Appendix A. Supplementary material

{supp}

# Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

# Acknowledgements

{ack}

# References

::: {{#refs}}
:::
"""
figs = s[s.index("# Figures\n"):]
widths = {1: "60%", 2: "60%", 3: "95%", 4: "60%", 5: "95%"}
fig_md = ["# Figures\n"]
for m in re.finditer(r"\*\*Figure (\d+)\. ([^*]+)\*\* ([^\n]+)", figs):
    n = int(m.group(1))
    fig_md.append(f"![**Figure {n}. {m.group(2)}** {m.group(3)}](figures/Figure{n}.pdf){{width={widths.get(n, '80%')}}}\n")
head = f"""---
title: "{title}"
---

**Ihor Kendiukhov**^a,\\*^

^a^ {aff}

^\\*^ Corresponding author. E-mail: {corr}

# Abstract

{abstract}

**Keywords:** {"; ".join(kw)}

"""
open(os.path.join(CBC, "manuscript_cbc.md"), "w").write(head + body + end + "\n" + "\n".join(fig_md))
open(os.path.join(CBC, "highlights.txt"), "w").write("Highlights\n\n" + "\n".join(f"• {h}" for h in hl) + "\n")
open(os.path.join(CBC, "title_abstract_keywords.txt"), "w").write(
    title + "\n" + re.sub(r"[*]", "", abstract).replace("^", "") + "\n" + "; ".join(kw) + "\n")
open(os.path.join(CBC, "credit_author_statement.txt"), "w").write(
    "CRediT authorship contribution statement\n\nIhor Kendiukhov: Conceptualization, Methodology, Software, Formal analysis, "
    "Investigation, Data curation, Visualization, Writing – original draft, Writing – review & editing.\n")
open(os.path.join(CBC, "declaration_of_interest.txt"), "w").write(
    "Declaration of interests\n\nThe author declares that he has no known competing financial interests or personal "
    "relationships that could have appeared to influence the work reported in this paper.\n\nIhor Kendiukhov\n")
open(os.path.join(CBC, "title_page.md"), "w").write(
    f"# {title}\n\n**Ihor Kendiukhov**^a,\\*^\n\n^a^ {aff}\n\n^\\*^ Corresponding author: Ihor Kendiukhov, {aff}. "
    f"E-mail: {corr}\n\n**Short title:** Contextual gene representations in single-cell models\n\n"
    f"**Keywords:** {'; '.join(kw)}\n\n**Word count of the abstract:** {n_words}\n")
print(f"manuscript_cbc.md written; abstract {n_words} words; {len(hl)} highlights (max {max(map(len, hl))} chars); "
      f"{len(kw)} keywords")
