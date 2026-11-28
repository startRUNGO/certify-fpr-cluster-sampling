# P5a figure scripts: input loaders for the public repository (patch 2026-11-28).
# load_master accepts the public redacted master results/results_master.json.gz (or its uncompressed copy) and the internal masters the scripts were rendered from.
#   The redaction changed only path strings in `sources` and `files_sha256` file names (REDACTIONS.tsv); every value a figure reads is unchanged.
# need_file stops with exit code 3 and a plain message when a figure needs a file that is not in this repository (available from the corresponding author on request).
import gzip, hashlib, json, sys
from pathlib import Path

PUBLIC_MASTER_SHA = "4b29688492be62d98512d220cd6c8a4cd3d11dad2138ad5c4cc9ab6b585ff8b5"
ORIGINAL_MASTER_SHA = "131b06d0d0aa8ebeaba6f2d8c3fb102494106ca0182fbdb4d4b40c7913547489"
ON_REQUEST = "This file is not in the public repository; it is available from the corresponding author on reasonable request."


def stop(msg):
    print("GATE (exit 3, nothing written):", msg, file=sys.stderr)
    raise SystemExit(3)


def load_master(path, internal_pins=()):
    """Returns (master dict, sha256 of the uncompressed JSON).  Accepts the public master and the listed internal pins."""
    p = Path(path)
    if not p.is_file(): stop(f"master not found: {p}")
    raw = p.read_bytes()
    if p.suffix == ".gz": raw = gzip.decompress(raw)
    h = hashlib.sha256(raw).hexdigest()
    if h not in (PUBLIC_MASTER_SHA, ORIGINAL_MASTER_SHA) + tuple(internal_pins): stop(f"master sha {h[:16]} is neither the public master nor a pinned internal master")
    return json.loads(raw), h


def need_file(path, want_sha, what):
    """Stops with a plain message if a non-public input is missing or differs from its pin."""
    if path is None or not Path(path).is_file(): stop(f"{what} is needed for this figure. {ON_REQUEST}")
    h = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if h != want_sha: stop(f"{what}: sha {h[:16]} differs from the pin {want_sha[:16]}")
    return Path(path)
