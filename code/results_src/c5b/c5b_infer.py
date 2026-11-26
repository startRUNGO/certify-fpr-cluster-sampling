"""C5b inference (PREREG_P5_C5B_FRAMING_AUG_20261026.md and errata 1-3). Written like the reviewed lh_infer.py (preallocated uint8 cache filled from pool.imap, PR-OD front end only, batch 512).
For each of the 10 C5b checkpoints (c5b-f_s41..45, c5b-fs_s41..45) it writes the two-class logits (float32, n x 2, same layout as the frozen A2 files) on:
  nih_validation  the NIH VALIDATION fold only (asserted; the locked test fold is never opened): source threshold (Youden J, rule of p5_eval_phase1_iter1.py) and the fold's AUROC go to <stem>_eval.json (a file; nothing is printed)
  bbbc041_all     all 85,486 BBBC041 crops in the order of iter1's manifest.csv (the order p5_eval_phase1_iter1.py used for ids_iter1.npz); names / y / group are compared with --ids-iter1
  mh              MP-IDB host-cell set v1 (7,086 crops, manifest sha 9c117d4f...), the layout of mh_infer.py
  lh              rebuilt Lacuna-B crops_v3_hc (52,859 crops, manifest sha b61c2671...), the layout of lh_infer.py
  lb_same         the rows of ids_lb_same.npz (--ids-lb-same) read from the OLD crops_v2 folder (1,517 positives + 51,342 negatives), in the row order of that file (mechanism reading)
Besides the logits each site gets ids_<site>.npz (the layout of the frozen A2 ids files). No metric of any development site is computed or printed. Writes only into a new directory (refuses to overwrite).
Usage: python3 c5b_infer.py --out <new dir> [--ckpt-dir ~/P5_StainRef/runs/c5b/ckpt] --ids-iter1 <.../iter1_20261021/infer/ids_iter1.npz> --ids-lb-same <ids_lb_same.npz> [--sites a,b] [--limit N] [--manifests-only]"""
import argparse, csv, hashlib, json, os, sys, time
from pathlib import Path
import numpy as np
from PIL import Image

HOME = Path.home()
sys.path.insert(0, os.environ.get("P5_CODE", str(HOME / "P5_StainRef" / "code")))
IMG = 128
SITES_ALL = ["nih_validation", "bbbc041_all", "mh", "lh", "lb_same"]
SEEDS = [41, 42, 43, 44, 45]; ARMS = ["c5b-f", "c5b-fs"]
GIEMSA = [[[[0.15]], [[-0.70]], [[0.70]]], [[[0.70]], [[0.15]], [[-0.70]]], [[[0.58]], [[0.58]], [[0.58]]]]
EXPECT = {"nih_validation": (4081,), "bbbc041_all": (85486, 2452, 83034, 1328), "mh": (7086, 553, 6533, 104), "lh": (52859, 1517, 51342, 483), "lb_same": (52859, 1517, 51342, 483)}
MH_SHA = "9c117d4f9c142e166b91c80b0ca5c1c59e9b3ce131a3bf4fbff2b54e1eda09c4"; LH_SHA = "b61c2671cc5423f124a10bc717d8607f59b7f8d6f179dde5c9ae410c3d535957"
ROOTS = dict(nih=os.environ.get("C5B_NIH_ROOT", str(HOME / "MedPrune" / "datasets" / "malaria" / "cell_images")),
             bbbc=os.environ.get("C5B_BBBC_ROOT", str(HOME / "P5_StainRef" / "iter1_20261021" / "data" / "bbbc041_allcells")),
             mh=os.environ.get("C5B_MH_ROOT", "<DATA>/p6_pathotrunk/mpidb_hostcells_v1"), lh=os.environ.get("C5B_LH_ROOT", "<DATA>/p5_e7_lacuna/crops_v3_hc"),
             old=os.environ.get("C5B_OLD_ROOT", "<DATA>/p5_e7_lacuna/crops_v2"))
SMOKE = bool(os.environ.get("C5B_SMOKE"))          # local synthetic smoke only: the real set sizes and manifest shas are not asserted


def sha_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def counts(y, field):
    return (len(y), int(np.sum(y)), int(np.sum(y == 0)), len(set(field)))


def site_rows(site, a):
    """returns dict(paths, y, field, ids=dict of arrays to store) for one site and asserts the set sizes and (where an archive exists) the identity with the frozen A2 ids."""
    if site == "nih_validation":
        from p5stage.data import scan_classification_dataset, split_cases_greedy
        samples = scan_classification_dataset(Path(ROOTS["nih"])); split_map = split_cases_greedy(samples, seed=42)
        items = [s for s in samples if split_map[s.case_id] == "validation"]; items.sort(key=lambda s: s.path)       # the order of nih_fold("validation") in p5_eval_phase1_iter1.py
        assert all(split_map[s.case_id] == "validation" for s in items), "a non-validation sample is in the NIH list"
        paths = [str(Path(ROOTS["nih"]) / s.path) for s in items]; y = np.array([s.label for s in items], dtype=int); field = np.array([s.case_id for s in items])
        names = np.array([Path(p).name for p in paths])
        assert SMOKE or len(y) == EXPECT[site][0], "NIH validation fold size"
        if a.ids_iter1 and not SMOKE:
            z = np.load(a.ids_iter1, allow_pickle=True)
            assert [str(x) for x in z["nih_validation_name"]] == [str(x) for x in names] and np.array_equal(z["nih_validation_y"], y), "NIH validation fold differs from ids_iter1"
        return dict(paths=paths, y=y, field=field, ids=dict(name=names, y=y, field=field))
    if site == "bbbc041_all":
        root = Path(ROOTS["bbbc"]); mf = list(csv.DictReader(open(root / "manifest.csv")))
        paths = [str(root / r["file"]) for r in mf]; y = np.array([int(r["label"]) for r in mf], dtype=int); field = np.array([r["image_uuid"] for r in mf])
        assert SMOKE or counts(y, field) == EXPECT[site], "BBBC-ALL set size assertion failed: %s" % (counts(y, field),)
        if a.ids_iter1 and not SMOKE:
            z = np.load(a.ids_iter1, allow_pickle=True)
            assert np.array_equal(z["bbbc041_all_y"], y) and [str(x) for x in z["bbbc041_all_group"]] == [str(x) for x in field] and [str(x) for x in z["bbbc041_all_name"]] == [Path(p).name for p in paths], "BBBC-ALL order/labels/groups differ from ids_iter1"
        return dict(paths=paths, y=y, field=field, ids=dict(crop_file=np.array([r["file"] for r in mf]), field=field, y=y))
    if site == "mh":
        root = Path(ROOTS["mh"])
        if not SMOKE: assert sha_file(root / "manifest.csv") == MH_SHA, "MP-IDB manifest sha mismatch"
        raw = list(csv.DictReader(open(root / "manifest.csv")))
        rows = [{"file": f"{r['class']}/{r['field']}_c{r['instance']}.png", "label": 1 if r["class"] == "Parasitized" else 0, "field_id": r["field"], "kind": "pos" if r["class"] == "Parasitized" else "neg"} for r in raw]
        paths = [str(root / r["file"]) for r in rows]; y = np.array([r["label"] for r in rows]); field = np.array([r["field_id"] for r in rows])
        assert SMOKE or counts(y, field) == EXPECT[site], "MP-IDB set size assertion failed: %s" % (counts(y, field),)
        return dict(paths=paths, y=y, field=field, ids=dict(crop_file=np.array([r["file"] for r in rows]), field=field, kind=np.array([r["kind"] for r in rows]), y=y))
    if site == "lh":
        root = Path(ROOTS["lh"])
        if not SMOKE: assert sha_file(root / "manifest.csv") == LH_SHA, "crops_v3_hc manifest sha mismatch"
        raw = list(csv.DictReader(open(root / "manifest.csv")))
        assert all(r["subsite"] == "B" and r["tuning"] == "0" and (r["kind"] == "pos") == (r["label"] == "1") for r in raw), "manifest rows are not all held-out Lacuna-B with consistent kind/label"
        paths = [str(root / r["file"]) for r in raw]; y = np.array([int(r["label"]) for r in raw]); field = np.array([r["field_id"] for r in raw])
        assert SMOKE or counts(y, field) == EXPECT[site], "rebuilt Lacuna-B set size assertion failed: %s" % (counts(y, field),)
        return dict(paths=paths, y=y, field=field, ids=dict(crop_file=np.array([r["file"] for r in raw]), field=field, kind=np.array([r["kind"] for r in raw]), y=y))
    if site == "lb_same":
        z = np.load(a.ids_lb_same, allow_pickle=True); root = Path(ROOTS["old"])
        paths = [str(root / str(f)) for f in z["crop_file"]]; y = np.array(z["y"], dtype=int); field = np.array(z["field"])
        assert SMOKE or counts(y, field) == EXPECT[site], "lb_same set size assertion failed: %s" % (counts(y, field),)
        return dict(paths=paths, y=y, field=field, ids={k: z[k] for k in z.files})
    raise SystemExit("unknown site " + site)


def load_one(args):
    p, pr = args; arr = np.asarray(Image.open(p).convert("RGB"), dtype=np.uint8)
    if pr:
        from p5stage.prod import legacy_float_pr_od
        arr = legacy_float_pr_od(arr)
    return np.asarray(Image.fromarray(arr).resize((IMG, IMG), Image.BILINEAR), dtype=np.uint8)


def cache(paths, pr=True):
    """one preallocated (N,3,IMG,IMG) uint8 array filled in order from pool.imap (as lh_infer.py); returns a numpy array."""
    from multiprocessing import Pool
    out = np.empty((len(paths), 3, IMG, IMG), dtype=np.uint8)
    with Pool(int(os.environ.get("NPROC", "4"))) as pool:
        for i, a in enumerate(pool.imap(load_one, [(p, pr) for p in paths], chunksize=64)): out[i] = a.transpose(2, 0, 1)
    return out


def youden_threshold(y, p):
    from sklearn.metrics import roc_curve
    fpr, tpr, th = roc_curve(y, p); return float(th[(tpr - fpr).argmax()])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--ckpt-dir", default=str(HOME / "P5_StainRef" / "runs" / "c5b" / "ckpt"))
    ap.add_argument("--ids-iter1", default=None); ap.add_argument("--ids-lb-same", default=None)
    ap.add_argument("--sites", default=",".join(SITES_ALL)); ap.add_argument("--limit", type=int, default=0, help="smoke test: first N rows per site")
    ap.add_argument("--manifests-only", action="store_true", help="build every site list and assert counts / identities; no torch, no images, nothing written")
    a = ap.parse_args(); sites = a.sites.split(",")
    assert set(sites) <= set(SITES_ALL)
    if "bbbc041_all" in sites or "nih_validation" in sites: assert a.ids_iter1, "--ids-iter1 is required for nih_validation / bbbc041_all"
    if "lb_same" in sites: assert a.ids_lb_same, "--ids-lb-same is required for lb_same"
    data = {s: site_rows(s, a) for s in sites}
    for s, d in data.items(): print(s, "n", len(d["y"]), "pos", int(d["y"].sum()), "fields", len(set(d["field"])), flush=True)
    if a.manifests_only: print("manifests ok"); return
    out = Path(a.out)
    if out.exists(): raise SystemExit("output dir exists: %s" % out)
    ck = sorted(Path(a.ckpt_dir).glob("*.pth")); stems = [p.stem for p in ck]
    assert stems == sorted("%s_s%d" % (g, s) for g in ARMS for s in SEEDS), "the checkpoint folder must hold exactly the 10 C5b checkpoints: %s" % stems
    out.mkdir(parents=True)
    import torch
    torch.set_num_threads(int(os.environ.get("P5_TORCH_THREADS", "0")) or torch.get_num_threads())
    import importlib.util
    spec = importlib.util.spec_from_file_location("malarianet", os.environ["P5_MODEL_PY"]); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    DEV = "cuda"; mean = torch.tensor([0.485, 0.456, 0.406]).view(1, 3, 1, 1).to(DEV); std = torch.tensor([0.229, 0.224, 0.225]).view(1, 3, 1, 1).to(DEV)

    def build_model(seed):
        torch.manual_seed(seed); m = mod.MalariaNet()
        with torch.no_grad(): m.scfe.stain_proj.weight.copy_(torch.tensor(GIEMSA, dtype=m.scfe.stain_proj.weight.dtype))
        m.scfe.stain_proj.weight.requires_grad = False; return m

    def strip(state):
        if not any(k.startswith("base.") for k in state): return state
        return {k[5:]: v for k, v in state.items() if k.startswith("base.")}

    @torch.no_grad()
    def infer(model, c, bs=512):
        lg = []
        for i in range(0, len(c), bs):
            x = torch.from_numpy(c[i:i + bs]).to(DEV, non_blocking=True).float().div_(255.0); x = (x - mean) / std; o = model(x); o = o[0] if isinstance(o, tuple) else o; lg.append(o.float().cpu())
        return torch.cat(lg).numpy()

    def get_model(pth):
        seed = int(pth.stem.rsplit("_s", 1)[1])
        blob = torch.load(pth, map_location="cpu", weights_only=False); m = build_model(seed); m.load_state_dict(strip(blob["model_state_dict"]), strict=True); return m.to(DEV).eval()

    meta = dict(script_sha256=sha_file(__file__), ckpt_sha256={p.stem: sha_file(p) for p in ck}, sites={}, torch=torch.__version__, cuda=torch.version.cuda, threads=torch.get_num_threads())
    t0 = time.time()
    for s in sites:
        d = data[s]; paths, y = d["paths"], d["y"]
        if a.limit: paths, y = paths[:a.limit], y[:a.limit]
        assert all(os.path.exists(p) for p in paths), "missing crop file in site " + s
        np.savez_compressed(out / ("ids_%s.npz" % s), **({k: (v[:a.limit] if a.limit else v) for k, v in d["ids"].items()}))
        c = cache(paths, True); print("cached", s, len(paths), "%.0fs" % (time.time() - t0), flush=True)
        for pth in ck:
            stem = pth.stem; model = get_model(pth); lg = infer(model, c)
            assert lg.shape == (len(paths), 2) and np.isfinite(lg).all(), "non-finite or mis-shaped logits"
            np.savez_compressed(out / ("%s_%s_logits.npz" % (stem, s)), logits=lg)
            if s == "nih_validation":
                from sklearn.metrics import roc_auc_score
                p = torch.softmax(torch.from_numpy(lg), 1)[:, 1].numpy()                    # the probabilities of p5_eval_phase1_iter1.py infer()
                json.dump({"ckpt": stem, "ckpt_sha256": meta["ckpt_sha256"][stem], "source_threshold": {"from": "nih_validation", "rule": "Youden J", "value": youden_threshold(y, p)},
                           "nih_validation": {"auroc": float(roc_auc_score(y, p)), "n": int(len(y))}}, open(out / ("%s_eval.json" % stem), "w"), indent=1)
            print("DONE", stem, s, "%.0fs" % (time.time() - t0), flush=True); del model; torch.cuda.empty_cache()
        meta["sites"][s] = dict(n=int(len(paths)), peak_cuda_alloc_MiB=torch.cuda.max_memory_allocated() / 2 ** 20); del c
    json.dump(meta, open(out / "c5b_infer_meta.json", "w"), indent=1)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
