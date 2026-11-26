"""C5b memory trial (not part of the evaluation; written like lh_memtrial2.py). On the first N rows (default 5,000) of the BBBC-ALL manifest and ONE checkpoint (default c5b-f_s41) it builds the PR-OD cache with c5b_infer.cache() (the code of the real run),
runs the 512-image batches through the model on the GPU, and samples the Pss of the whole process tree (this process + pool workers, /proc/<pid>/smaps_rollup, every 0.2 s) and the torch GPU peak.
It then extrapolates linearly to the 85,486 rows of BBBC-ALL, the largest site (extrapolated = base + (peak - base) x 85486 / N, base = tree Pss after the imports), and prints the arithmetic size of the full cache (85,486 x 3 x 128 x 128 bytes).
No logits are saved, no score, label or metric is computed or printed; only memory figures. Usage: NPROC=4 python3 c5b_memtrial.py <ckpt path> [N=5000]  (env: P5_MODEL_PY, C5B_BBBC_ROOT, P5_TORCH_THREADS)"""
import os, sys, threading, time, resource, platform
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import c5b_infer as M
FULL = 85486; N = int(sys.argv[2]) if len(sys.argv) > 2 else 5000; G = 2 ** 30


def descendants(root):
    kids = {}
    for d in os.listdir("/proc"):
        if d.isdigit():
            try: kids.setdefault(int(open("/proc/%s/stat" % d).read().rsplit(")", 1)[1].split()[1]), []).append(int(d))
            except Exception: pass
    out, todo = [], [root]
    while todo:
        p = todo.pop(); out.append(p); todo += kids.get(p, [])
    return out


def tree_pss_gib():
    if not os.path.exists("/proc/self/smaps_rollup"): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if platform.system() == "Darwin" else 1024) / G     # local smoke only
    s = 0
    for p in descendants(os.getpid()):
        try:
            for l in open("/proc/%d/smaps_rollup" % p):
                if l.startswith("Pss:"): s += int(l.split()[1]); break
        except Exception: pass
    return s * 1024 / G


peak = [0.0]; stop = threading.Event()
def sampler():
    while not stop.is_set(): peak[0] = max(peak[0], tree_pss_gib()); time.sleep(0.2)


if __name__ == "__main__":
    ckpt = Path(sys.argv[1]); import csv
    root = Path(M.ROOTS["bbbc"]); mf = list(csv.DictReader(open(root / "manifest.csv")))[:N]; paths = [str(root / r["file"]) for r in mf]
    assert all(os.path.exists(p) for p in paths)
    import torch
    torch.set_num_threads(int(os.environ.get("P5_TORCH_THREADS", "0")) or torch.get_num_threads())
    base = tree_pss_gib(); t = threading.Thread(target=sampler, daemon=True); t.start()
    print("%d of %d rows; base (tree Pss after imports) %.2f GiB; NPROC %s" % (len(paths), FULL, base, os.environ.get("NPROC", "4")), flush=True)
    c = M.cache(paths, True); after_cache = peak[0]
    print("after the PR-OD cache: tree Pss peak so far %.2f GiB; cache array %.3f GiB (full set: %.2f GiB)" % (after_cache, c.nbytes / G, FULL * 3 * M.IMG * M.IMG / G), flush=True)
    import importlib.util
    spec = importlib.util.spec_from_file_location("malarianet", os.environ["P5_MODEL_PY"]); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    seed = int(ckpt.stem.rsplit("_s", 1)[1]); torch.manual_seed(seed); m = mod.MalariaNet()
    blob = torch.load(ckpt, map_location="cpu", weights_only=False); st = blob["model_state_dict"]; st = {k[5:]: v for k, v in st.items() if k.startswith("base.")} if any(k.startswith("base.") for k in st) else st
    m.load_state_dict(st, strict=True); m = m.to("cuda").eval(); mean = torch.tensor([0.485, 0.456, 0.406]).view(1, 3, 1, 1).cuda(); std = torch.tensor([0.229, 0.224, 0.225]).view(1, 3, 1, 1).cuda()
    with torch.no_grad():
        for i in range(0, len(c), 512): x = torch.from_numpy(c[i:i + 512]).cuda().float().div_(255.0); m((x - mean) / std)
    time.sleep(0.5); stop.set(); t.join()
    gpu = torch.cuda.max_memory_allocated() / 2 ** 20; scale = FULL / len(paths)
    print("after the inference: tree Pss peak %.2f GiB; torch GPU peak allocated %.0f MiB (batch-bound, does not scale with N)" % (peak[0], gpu), flush=True)
    est = base + (peak[0] - base) * scale
    print("extrapolation to %d rows: tree Pss = %.1f GiB (base %.2f + (peak %.2f - base) x %.2f); + context (CUDA) is inside the measured peak" % (FULL, est, base, peak[0], scale), flush=True)
    print("ALL DONE", flush=True)
