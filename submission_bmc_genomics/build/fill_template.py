"""Render a template (e.g. manuscript_bmc.template.md) by replacing every {{key}} with build/numbers.json[key].
Fails (exit 1) if any placeholder has no value. Usage: fill_template.py template.md out.md"""
import sys, os, re, json
HERE = os.path.dirname(os.path.abspath(__file__))
N = json.load(open(os.path.join(HERE, "numbers.json")))
src, dst = sys.argv[1], sys.argv[2]
s = open(src).read()
missing = sorted({k for k in re.findall(r"\{\{([A-Za-z0-9_]+)\}\}", s) if k not in N})
if missing:
    sys.exit(f"unresolved placeholders in {src}: {missing}")
s = re.sub(r"\{\{([A-Za-z0-9_]+)\}\}", lambda m: str(N[m.group(1)]), s)
open(dst, "w").write(s)
print(f"rendered {dst}")
