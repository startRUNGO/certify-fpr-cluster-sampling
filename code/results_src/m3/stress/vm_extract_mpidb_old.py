# P5 stress test (a), VM side: extract ONLY the key "mpidb_full_logits" from the iter1 A2 score files of seeds 42-45 (old MP-IDB, a development-site version).
# WRITTEN, NOT RUN, until 实验管理 (VM manager) approves the exact text.  Read-only on the sources; writes < 1 MB to a NEW UTC directory under the VM home.
# Rules (实验管理 2026-10-30):
#   - each source file's sha256 must equal the source_sha256 of a2_link/MANIFEST.tsv (pinned below, copied by program from MANIFEST); otherwise exit, nothing written;
#   - np.load is lazy per key: only z["mpidb_full_logits"] is accessed; no z.files, no dict(z), no iteration, no other key (the files also hold NIH keys);
#   - output: 4 npz (key "logits") + SHA256SUMS, opened with "x"; if the output directory exists, exit;
#   - prints only the output paths and sha256 lines.
# Run as: nice -n 19 python -I vm_extract_mpidb_old.py   (single process; does not touch <DATA>/lvzhonghao or <VM_HOME>/jianbo.huang)
import hashlib, os, sys, time
from pathlib import Path
import numpy as np

SRC = Path("<VM_HOME>/P5_StainRef/runs/iter1_20261021/infer")
PIN = {42: "2afcc4c3a0cc9ef4d846395405e9eb23b3bc85ba8f89e5100551edb924ee9d35", 43: "6a3461e068fef2273e7e79cf489ec075cf186749877439d02e4b9c780ddac0a6", 44: "f74cb101b1d7a9e84928d423f287b8c871b85044f37f31a62e9fc603a607086a", 45: "95ed320994c83907df68a3c98b486e81dabd047116739665762b4cf78048e83e"}
OUT = Path.home() / ("p5_mpidb_old_extract_" + time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()))


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def main():
    for sd, want in PIN.items():
        if sha(SRC / f"a2_s{sd}_scores.npz") != want: sys.exit(f"source sha differs for seed {sd}; nothing written")
    os.mkdir(OUT)                                                            # raises if it exists
    lines = []
    for sd in PIN:
        z = np.load(SRC / f"a2_s{sd}_scores.npz")
        arr = np.asarray(z["mpidb_full_logits"])                            # the only key accessed
        out = OUT / f"a2_s{sd}_mpidb_full_logits.npz"
        with open(out, "xb") as fh: np.savez(fh, logits=arr)
        lines.append(f"{sha(out)}  ./{out.name}")
    with open(OUT / "SHA256SUMS", "x") as fh: fh.write("\n".join(lines) + "\n")
    print(OUT); print("\n".join(lines))


if __name__ == "__main__":
    main()
