"""Add the title and the cover-letter findings to build/numbers.json (they are composed from the numbers), and after
the manuscript is rendered write title_abstract_keywords.txt. Usage: compose_extra.py numbers | form"""
import re, json, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
if sys.argv[1] == "numbers":
    N = json.load(open(os.path.join(HERE, "numbers.json")))
    t = open(os.path.join(ROOT, "manuscript_bmc.template.md")).read()
    N["TITLE"] = re.search(r'title: "(.*)"', t).group(1)
    N["COVER_FINDINGS"] = f"""The main results are:

1. **The contextual representation is real.** In MaxToki-217M, a gene's representation changes with cell type in a gene-specific way that reproduces across independent cells (cosine {N['excess_L4_r2']} against {N['diff_L4_r3']} for a shuffled control), is strongest in the first two layers, is not explained by the gene's rank in the cell, and replicates on an independent bone-marrow dataset. Training strengthens it (an untrained model of the same architecture keeps about {N['cm_rand_excess_pct']}% of the effect but only about a third of its functional organisation), and it is weakly aligned with the model's objective: genes that change more with context gain more from real context in next-gene prediction. All four models tested (MaxToki-217M and -1B, scGPT, STATE) contextualise genes strongly.
2. **Nothing context-specific beyond co-expression.** The context changes follow functional axes, but gene sets matched on size and co-expression give similar organisation (p = {N['pm_range']}). Two direct tests for context-specific gene function, run in MaxToki-217M, are negative; the only signal beyond co-expression is a small, static closeness of transcription factors to their curated targets. Steering shows that the model reads one functional direction (nuclear versus surface genes); a second direction has no consistent effect and a third acts in reverse. Reading a gene's changed role in a new cell type is therefore not reachable with the model tested.
3. **A model metric.** The functional organisation of contextual change differs sharply between architectures (z = {N['cm_maxtoki_fz_range']} for MaxToki, {N['cm_scgpt_fz_r1']} for scGPT, {N['cm_state_fz_r1']} for STATE), while all models contextualise genes strongly."""
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
