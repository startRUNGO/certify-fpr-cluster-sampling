# Local smoke for vm_extract_mpidb_old.py with FAKE source files (no real data).  WRITTEN, NOT RUN, until 实验管理 approves.
# (1) happy path: fake a2_s{42..45}_scores.npz holding mpidb_full_logits plus NIH-named keys; outputs are 4 npz with only "logits" equal to the source array, and a
#     SHA256SUMS whose lines match the files; the ONLY key accessed on any source is "mpidb_full_logits" (np.load wrapped to log key access; .files/iteration logged).
# (2) a source whose sha differs from its pin -> exit, output directory not created.
# (3) output directory already exists -> refused (os.mkdir raises), nothing overwritten.
# Usage: python -I -B vm_extract_mpidb_old_smoke_local.py <empty scratch dir>
import hashlib, importlib.util, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; T = Path(sys.argv[1]); assert T.is_dir() and not any(T.iterdir())
spec = importlib.util.spec_from_file_location("vx", HERE / "vm_extract_mpidb_old.py"); V = importlib.util.module_from_spec(spec); spec.loader.exec_module(V)
ok = True


def report(tag, cond, extra=""):
    global ok
    ok &= bool(cond); print(tag, "ok" if cond else "FAIL", extra, flush=True)


def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


src = T / "src"; src.mkdir(); rng = np.random.default_rng(0); truth = {}
for sd in (42, 43, 44, 45):
    a = rng.normal(size=(30, 2)); truth[sd] = a
    np.savez(src / f"a2_s{sd}_scores.npz", mpidb_full_logits=a, nih_locked_test_logits=np.ones((5, 2)), bbbc041_all_logits=np.zeros((7, 2)))
V.SRC = src; V.PIN = {sd: h(src / f"a2_s{sd}_scores.npz") for sd in (42, 43, 44, 45)}
ACC = []; _load = np.load


class Proxy:
    def __init__(self, z): self.z = z
    def __getitem__(self, k): ACC.append(k); return self.z[k]
    def __getattr__(self, k): ACC.append("ATTR:" + k); return getattr(self.z, k)
    def __iter__(self): ACC.append("ITER"); return iter(self.z)


V.np.load = lambda p, *a, **k: Proxy(_load(p, *a, **k))
V.OUT = T / "out1"
try: V.main(); e1 = None
except SystemExit as e: e1 = repr(e)
V.np.load = _load
sums = (T / "out1" / "SHA256SUMS").read_text().split("\n") if e1 is None else []
same = e1 is None and all(np.array_equal(np.load(T / "out1" / f"a2_s{sd}_mpidb_full_logits.npz")["logits"], truth[sd]) and
                          np.load(T / "out1" / f"a2_s{sd}_mpidb_full_logits.npz").files == ["logits"] for sd in truth)
sums_ok = e1 is None and all(f"{h(T / 'out1' / f'a2_s{sd}_mpidb_full_logits.npz')}  ./a2_s{sd}_mpidb_full_logits.npz" in sums for sd in truth)
report("(1) happy path, only one key accessed", same and sums_ok and set(ACC) == {"mpidb_full_logits"} and len(ACC) == 4, sorted(set(ACC)))
V.PIN[44] = "0"*64; V.OUT = T / "out2"
try: V.main(); e2 = None
except SystemExit as e: e2 = repr(e)
report("(2) sha mismatch refused", e2 is not None and "44" in e2 and not (T / "out2").exists(), e2)
V.PIN[44] = h(src / "a2_s44_scores.npz"); V.OUT = T / "out1"
try: V.main(); e3 = None
except FileExistsError as e: e3 = repr(e)
report("(3) existing output dir refused", e3 is not None, e3)
print("ALL PASS" if ok else "FAIL")
