#!/bin/bash
# Build the BMC Genomics manuscript: DOCX (main submission file) and a PDF preview, from manuscript_bmc.md.
# DOCX: Times New Roman 12 pt, double spacing, continuous line numbers, page numbers, numbered (Vancouver) refs.
set -e; cd "$(dirname "$0")/.."
export PATH="/opt/homebrew/bin:/Library/TeX/texbin:$PATH"
PANDOC_REFS="--citeproc --bibliography=references_bmc.bib --csl=build/biomed-central.csl"
# 1) every number from the result files -> build/numbers.json; 2) render the templates (fails on any unresolved key)
python3 build/make_numbers.py && python3 build/compose_extra.py numbers
python3 build/make_additional_file.py      # Additional file 1 (Excel) from the same result files
python3 build/fill_template.py manuscript_bmc.template.md manuscript_bmc.md
python3 build/fill_template.py cover_letter_bmc.template.md cover_letter_bmc.md
python3 build/compose_extra.py form
# 3) figures (no titles inside the graphic, BMC widths) -> figures/Figure1..5.pdf
(cd .. && FIGDIR=bmc_submission/figures FIG_TITLES=0 python3 ctx_figures.py >/dev/null)
for n in 1 2 3 4 5; do mv -f figures/ctx_fig$n.pdf figures/Figure$n.pdf; done
pandoc manuscript_bmc.md -f markdown -t docx $PANDOC_REFS --reference-doc=build/reference_bmc.docx -o build/_raw.docx
python3 build/postprocess_docx.py build/_raw.docx manuscript_bmc.docx && rm build/_raw.docx
# PDF preview (same text; double spacing, line numbers) for checking only -- the DOCX is the file to upload
pandoc manuscript_bmc.md -f markdown -t latex --standalone $PANDOC_REFS \
  -V documentclass=article -V fontsize=12pt -V geometry:"a4paper,margin=2.5cm" -V linestretch=2 \
  -V mainfont="Times New Roman" -V colorlinks=true -V urlcolor=NavyBlue -V linkcolor=NavyBlue -V citecolor=NavyBlue \
  -H build/lineno_header.tex -o build/manuscript_bmc_preview.tex
(cd build && latexmk -xelatex -interaction=nonstopmode manuscript_bmc_preview.tex >/dev/null 2>&1 || true; latexmk -c manuscript_bmc_preview.tex >/dev/null 2>&1 || true)
mv build/manuscript_bmc_preview.pdf manuscript_bmc_preview.pdf
echo "built manuscript_bmc.docx and manuscript_bmc_preview.pdf"
# cover letter (single spacing)
pandoc cover_letter_bmc.md -f markdown -t docx --reference-doc=build/reference_letter.docx -o cover_letter_bmc.docx
pandoc cover_letter_bmc.md -f markdown -t latex --standalone -V fontsize=11pt -V geometry:"a4paper,margin=2.5cm" \
  -V mainfont="Times New Roman" -V colorlinks=true -V urlcolor=NavyBlue -V pagestyle=empty -o build/cover_letter_bmc.tex
(cd build && latexmk -xelatex -interaction=nonstopmode cover_letter_bmc.tex >/dev/null 2>&1 || true; latexmk -c cover_letter_bmc.tex >/dev/null 2>&1 || true)
mv build/cover_letter_bmc.pdf cover_letter_bmc.pdf
echo "built cover_letter_bmc.docx and cover_letter_bmc.pdf"
