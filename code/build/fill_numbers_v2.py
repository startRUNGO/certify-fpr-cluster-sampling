"""Replace [[id :: fmt]] placeholders in sections_src/*.tex with values from results_master.json -> sections/*.tex.
fmt: d (integer; four digits plain, five or more digits with a LaTeX thin space, 37\\,899, in text and in math), f2/f3/f4 (fixed decimals), pct1 (value*100, 1 decimal, no % sign), raw.
Unknown id => error. Also lists every number used with its source in used_numbers.tsv."""
import json, re, sys
from pathlib import Path
P = Path(__file__).resolve().parents[1]; M = json.load(open(P/"results_master.json")); V, S = M["values"], M["sources"]
pat = re.compile(r"\[\[\s*(.+?)\s*::\s*(\w+)\s*\]\]"); used = {}; bad = []
miss = re.compile(r"\[\[\s*MISSING:\s*(.+?)\s*\]\]")   # author-supplied items: one visible marker in the PDF, never silently dropped
def fmt(v, f, inmath=False):
    if f == "d":
        n = int(round(float(v))); return str(n) if abs(n) < 10000 else f"{n:,}".replace(",", "\\,")
    if f == "sf3": s = f"{float(v):.3f}"; return "+" + s if float(s) > 0 else s   # signed: positive values get "+", negatives keep "-" (turned into $-$ by sub())
    if f.startswith("f"): return f"{float(v):.{int(f[1:])}f}"
    if f == "pct1": return f"{100*float(v):.1f}"
    if f == "sci":   # scientific notation as math, e.g. 1e-06 -> $10^{-6}$, 3e-05 -> $3\times10^{-5}$
        m_, e_ = f"{float(v):.0e}".split("e"); b_ = f"10^{{{int(e_)}}}" if m_ == "1" else f"{m_}\\times10^{{{int(e_)}}}"; return b_ if inmath else "$" + b_ + "$"
    if f == "raw": return str(v)
    raise ValueError(f)
def sub(m):
    i, f = m.group(1), m.group(2)
    inmath = m.string[:m.start()].replace("\\$", "").count("$") % 2 == 1
    if i not in V: bad.append(i); return "??"
    used[i] = (V[i], S[i]); s = fmt(V[i], f, inmath)
    return s if inmath or not s.startswith("-") else "$-$" + s[1:]
for src in sorted((P/"sections_src").glob("*.tex")):
    t = miss.sub(lambda m: r"\textbf{[MISSING: " + m.group(1) + "]}", src.read_text()); t = pat.sub(sub, t); (P/"sections"/src.name).write_text(t)
if bad: print("UNKNOWN IDS:", sorted(set(bad))); sys.exit(1)
with open(P/"used_numbers.tsv", "w") as f:
    f.write("id\tvalue\tsource_file\tpath_or_derivation\n")
    # v2 (assemble draft 9 r3): derived keys may name several files ('files', sim50/platt_fallback) or none ('from_keys', sim50_derived)
    for i, (v, s) in sorted(used.items()):
        src_ = s["file"] if "file" in s else ";".join(s["files"]) if "files" in s else "results_master keys " + ",".join(s["from_keys"][:3]) + (",..." if len(s["from_keys"]) > 3 else "")
        f.write(f"{i}\t{v}\t{src_}\t{'|'.join(s['path']) if 'path' in s else s['derived']}\n")
print("filled", len(used), "distinct numbers")
