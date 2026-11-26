"""E8 PlasmoCount inference per PREREG_E8_PLASMOCOUNT_20261023.md (sha256 2ce88dfd...). Same pipeline as e7_infer.py (sha256 9c666a9f...): 35 checkpoints, batch 512."""
import os, sys, csv, json, time, numpy as np, torch
from pathlib import Path
from multiprocessing import Pool
from PIL import Image
HOME = Path.home(); sys.path.insert(0, str(HOME/"P5_StainRef"/"code"))
from p5stage.prod import legacy_float_pr_od
import importlib.util
spec = importlib.util.spec_from_file_location("malarianet", os.environ["P5_MODEL_PY"]); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
GIEMSA = [[[[0.15]], [[-0.70]], [[0.70]]], [[[0.70]], [[0.15]], [[-0.70]]], [[[0.58]], [[0.58]], [[0.58]]]]
IMG = 128; DEV = "cuda"; OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True); CROPS = "<DATA>/p5_e8_candidate/crops"
IMNET_MEAN = torch.tensor([0.485, 0.456, 0.406]).view(1, 3, 1, 1); IMNET_STD = torch.tensor([0.229, 0.224, 0.225]).view(1, 3, 1, 1)
def build_model(seed):
    torch.manual_seed(seed); m = mod.MalariaNet()
    with torch.no_grad(): m.scfe.stain_proj.weight.copy_(torch.tensor(GIEMSA, dtype=m.scfe.stain_proj.weight.dtype))
    m.scfe.stain_proj.weight.requires_grad = False; return m
def strip(state):
    if not any(k.startswith("base.") for k in state): return state
    return {k[5:]: v for k, v in state.items() if k.startswith("base.")}
def load_one(args):
    p, pr = args; arr = np.asarray(Image.open(p).convert("RGB"), dtype=np.uint8)
    if pr: arr = legacy_float_pr_od(arr)
    return np.asarray(Image.fromarray(arr).resize((IMG, IMG), Image.BILINEAR), dtype=np.uint8)
def cache(paths, pr):
    with Pool(int(os.environ.get("NPROC", "4"))) as pool: res = pool.map(load_one, [(p, pr) for p in paths], chunksize=64)
    return torch.from_numpy(np.stack(res)).permute(0, 3, 1, 2).contiguous()
@torch.no_grad()
def infer(model, c, bs=512):
    lg = []; mean = IMNET_MEAN.to(DEV); std = IMNET_STD.to(DEV)
    for i in range(0, len(c), bs):
        x = c[i:i+bs].to(DEV, non_blocking=True).float().div_(255.0); x = (x-mean)/std; o = model(x); o = o[0] if isinstance(o, tuple) else o; lg.append(o.float().cpu())
    return torch.cat(lg).numpy()
if __name__ == "__main__":
    rows = list(csv.DictReader(open(CROPS+"/manifest.csv"))); paths = [f"{CROPS}/{r['file']}" for r in rows]; y = np.array([int(r["label"]) for r in rows]); print("crops", len(rows), "pos", int(y.sum()), "neg", int((y == 0).sum()), flush=True)
    np.savez_compressed(OUT/"ids_plasmo.npz", crop_file=np.array([r["file"] for r in rows]), field=np.array([r["field_id"] for r in rows]), kind=np.array([r["kind"] for r in rows]), y=y)
    t0 = time.time(); caches = {False: cache(paths, False), True: cache(paths, True)}; print("cached", time.time()-t0, flush=True)
    for pth in sorted((HOME/"P5_StainRef/runs/main/ckpt").glob("*.pth")):
        stem = pth.stem; arm = stem.rsplit("_s", 1)[0]; seed = int(stem.rsplit("_s", 1)[1]); pr = arm in ("a2", "a3", "a4")
        blob = torch.load(pth, map_location="cpu", weights_only=False); model = build_model(seed); model.load_state_dict(strip(blob["model_state_dict"]), strict=True); model = model.to(DEV).eval()
        lg = infer(model, caches[pr]); np.savez_compressed(OUT/f"{stem}_plasmo_logits.npz", logits=lg)
        from sklearn.metrics import roc_auc_score
        auc = float(roc_auc_score(y, lg[:, 1]-lg[:, 0])); json.dump({"ckpt": stem, "auroc": auc, "peak_cuda_alloc_MiB": torch.cuda.max_memory_allocated()/2**20, "n": int(len(y))}, open(OUT/f"{stem}_plasmo_meta.json", "w"))
        print(f"DONE {stem} auroc {auc:.4f} ({time.time()-t0:.0f}s)", flush=True); del model; torch.cuda.empty_cache()
    print("ALL DONE", flush=True)
