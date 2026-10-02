"""Add the title and the cover-letter findings to build/numbers.json (they are composed from the numbers), and after
the manuscript is rendered write title_abstract_keywords.txt. Usage: compose_extra.py numbers | form"""
import re, json, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
if sys.argv[1] == "numbers":
    N = json.load(open(os.path.join(HERE, "numbers.json")))
    t = open(os.path.join(ROOT, "manuscript_bmc.template.md")).read()
    N["TITLE"] = re.search(r'title: "(.*)"', t).group(1)
    N["COVER_FINDINGS"] = f"""The main results are:

1. **The contextual representation is real, in every model tested, but not unique to models.** MaxToki-217M, MaxToki-1B, scGPT and STATE, run on identical cells, all change a gene's representation with cell type in a gene-specific way that reproduces across independent cells (agreement between independent halves of the cells {N['mm_excess_range']} above a shuffled control) and is not explained by the gene's rank in the cell. In MaxToki-217M it replicates on an independent bone-marrow dataset. Untrained models keep {N['cm_randpct_range']}% of it, and a representation built from expression alone, with no model, also shows it ({N['mm_exprC_excess']} when each gene vector averages 50 cells, as in the models).
2. **What the models add is organisation by co-expression.** In MaxToki the change follows functional axes (for example nuclear versus surface genes); STATE shows weaker organisation on that axis and scGPT little. In no model is this organisation significantly stronger than that of gene sets matched on size and co-expression. In every model the change follows co-expressed gene sets; in three representations built from expression alone, on the same cells, it does not. {N['COVER_STEER']} In both MaxToki sizes, genes that change more gain slightly more from their own context in next-gene prediction.
3. **No context-specific gene function beyond co-expression.** Two direct tests, run in all four models, are negative, although they could detect only clear effects. The only signal beyond co-expression is a static closeness of transcription factors to their curated targets that does not change with cell type; in STATE it is already present in the untrained model, which shares the trained model's protein-sequence gene table. The functional organisation of contextual change separates some architectures but does not yet rank them reliably."""
    json.dump(N, open(os.path.join(HERE, "numbers.json"), "w"), indent=1, ensure_ascii=False)
else:
    s = open(os.path.join(ROOT, "manuscript_bmc.md")).read()
    title = re.search(r'title: "(.*)"', s).group(1)
    ab = s[s.index("# Abstract") + len("# Abstract"):s.index("**Keywords:**")].strip()
    ab = re.sub(r"\^(.*?)\^", lambda m: m.group(1), ab).replace(" × 10−", " × 10^-")
    paras = [re.sub(r"\*\*(.*?)\*\*", r"\1", p.strip()) for p in ab.split("\n\n") if p.strip()]
    kw = s[s.index("**Keywords:**") + len("**Keywords:**"):].split("\n")[0].strip()
    open(os.path.join(ROOT, "title_abstract_keywords.txt"), "w").write(title + "\n" + " ".join(paras) + "\n" + kw + "\n")
    words = len(re.sub(r"\*\*|\^", "", s[s.index("# Abstract") + 10:s.index("**Keywords:**")]).split())
    print(f"abstract: {words} words incl. section labels (limit 350)")
    if words > 350:
        sys.exit("abstract over 350 words")
