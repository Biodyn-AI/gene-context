"""Write the general (preprint-style) manuscript ../PAPER_context_representation.md from the rendered BMC manuscript, so
both versions share one source. Differences: plain title/author header for make_paper_pdf.sh, figures embedded where their
legends are (from ../figures/ctx_figN.pdf), no 'Article type' line."""
import re, os
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); RG = os.path.dirname(ROOT)
s = open(os.path.join(ROOT, "manuscript_bmc.md")).read()
title = re.search(r'title: "(.*)"', s).group(1)
body = s[s.index("# Abstract"):]
aff = re.search(r"\^1\^ (.*)\n", s).group(1)
corr = re.search(r"Correspondence: (\S+)", s).group(1)
head = f"# {title}\n\n**Ihor Kendiukhov · 2026**\n\n{aff}. Correspondence: {corr}\n\n"
# embed figures: replace each legend paragraph in the '# Figures' section by an image with that caption
i = body.index("# Figures\n"); figs = body[i:]; body = body[:i]
widths = {1: "60%", 2: "60%", 3: "95%", 4: "60%", 5: "95%"}
out = ["# Figures\n"]
for m in re.finditer(r"\*\*Figure (\d+)\. ([^*]+)\*\* ([^\n]+)", figs):
    n = int(m.group(1))
    out.append(f"![**Figure {n}. {m.group(2)}** {m.group(3)}](figures/ctx_fig{n}.pdf){{width={widths[n]}}}\n")
open(os.path.join(RG, "PAPER_context_representation.md"), "w").write(head + body + "\n".join(out))
print("wrote PAPER_context_representation.md")
