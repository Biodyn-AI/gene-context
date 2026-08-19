#!/bin/bash
set -e; cd "$(dirname "$0")"
export PATH="/opt/homebrew/bin:/Library/TeX/texbin:$PATH"
MD="${1:-PAPER_context_representation.md}"; BASE="${MD%.md}"
python3 - "$MD" > /tmp/_ppm.txt <<'PY'
import sys,re
raw=open(sys.argv[1],encoding="utf-8").read()
raw="".join(c for c in raw if c in "\n\t" or ord(c)>=32)
lines=raw.split("\n"); title=lines[0].lstrip("# ").strip(); author=""; bs=1
for i in range(1,6):
    if re.fullmatch(r"\*\*.*\*\*",lines[i].strip()): author=lines[i].strip().strip("*").strip(); bs=i+1; break
open("/tmp/_ppb.md","w",encoding="utf-8").write("\n".join(lines[bs:]).lstrip("\n")); print(title); print(author)
PY
TITLE=$(sed -n '1p' /tmp/_ppm.txt); AUTHOR=$(sed -n '2p' /tmp/_ppm.txt)
pandoc /tmp/_ppb.md -f markdown -t latex --standalone --citeproc --bibliography=paper_ctx.bib \
  -M title="$TITLE" -M author="$AUTHOR" -M date="2026" -M reference-section-title="References" \
  -V documentclass=article -V fontsize=11pt -V geometry:"a4paper,margin=2.2cm" \
  -V colorlinks=true -V linkcolor=NavyBlue -V citecolor=NavyBlue -V urlcolor=NavyBlue \
  -o "$BASE.tex"
latexmk -xelatex -interaction=nonstopmode "$BASE.tex" >/dev/null 2>&1 || xelatex -interaction=nonstopmode "$BASE.tex" >/dev/null 2>&1
latexmk -c "$BASE.tex" >/dev/null 2>&1 || true
echo "wrote $BASE.pdf ($(grep -oE '[0-9]+ pages' "$BASE.log" 2>/dev/null | tail -1))"
