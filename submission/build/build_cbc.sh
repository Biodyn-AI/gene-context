#!/bin/bash
# Build the Computational Biology and Chemistry package from the rendered BMC manuscript (run build_bmc.sh first).
# Review PDF: single column, 1.5 spacing, continuous line numbers, figures with captions, numbered references.
set -e; cd "$(dirname "$0")/.."
export PATH="/opt/homebrew/bin:/Library/TeX/texbin:$PATH"
BMC=../bmc_submission
python3 build/make_cbc.py
mkdir -p figures && cp $BMC/figures/Figure*.pdf figures/
REFS="--citeproc --bibliography=$BMC/references_bmc.bib --csl=$BMC/build/biomed-central.csl"
pandoc manuscript_cbc.md -f markdown -t latex --standalone $REFS \
  -V documentclass=article -V fontsize=11pt -V geometry:"a4paper,margin=2.5cm" -V linestretch=1.5 \
  -V mainfont="Times New Roman" -V colorlinks=true -V urlcolor=NavyBlue -V linkcolor=NavyBlue -V citecolor=NavyBlue \
  -H $BMC/build/lineno_header.tex -o manuscript.tex
latexmk -xelatex -interaction=nonstopmode manuscript.tex >/dev/null 2>&1 || true
latexmk -c manuscript.tex >/dev/null 2>&1 || true
test -s manuscript.pdf || { echo "manuscript.pdf not built"; exit 1; }
# title page and cover letter
python3 $BMC/build/fill_template.py cover_letter_cbc.template.md cover_letter.md
for doc in title_page cover_letter; do
  pandoc $doc.md -f markdown -t docx --reference-doc=$BMC/build/reference_letter.docx -o $doc.docx
  pandoc $doc.md -f markdown -t latex --standalone -V fontsize=11pt -V geometry:"a4paper,margin=2.5cm" \
    -V mainfont="Times New Roman" -V colorlinks=true -V urlcolor=NavyBlue -V pagestyle=empty -o build/$doc.tex
  (cd build && latexmk -xelatex -interaction=nonstopmode $doc.tex >/dev/null 2>&1 || true; latexmk -c $doc.tex >/dev/null 2>&1 || true)
  mv build/$doc.pdf $doc.pdf
done
rm -f build/*.tex
echo "built the CBC package: manuscript.pdf/.tex, title_page, cover_letter, highlights, statements, figures"
