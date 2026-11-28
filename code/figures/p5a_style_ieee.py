# P5a shared figure style, IEEE Transactions version (p5a_style v2.0).  Style only: no data, no keys, no gates on data, no file writes.
# v2.0 (2026-11-24): IEEE style.
#   Kept from P4: serif Times New Roman with stix mathtext; 8 pt text, 7.5 pt ticks, 7 pt legend (nothing under 7 pt); axes 0.6; ticks inward; no top / right spines;
#   no grid (y-grid only on request, alpha 0.16); Okabe-Ito; colour + marker for every series; HERO / MUTED emphasis; pdf.fonttype 42; 600 dpi; constrained_layout.
#   Changed from P4: widths are the BJ widths (full 139 mm, half 67 mm), not IEEE 88.9 / 181.9 mm; HERO colours are the P5a role colours (S3 orange, SPLIT6 blue,
#   Z+ green, limits vermillion), not P4's HERO blue / orange; savefig.bbox is NOT "tight" (a tight crop would make the PDF narrower than 139 / 67 mm and LaTeX
#   would rescale it); render() checks the PDF /MediaBox width instead.
# The v1.1 module p5a_style.py stays unchanged (the approved v1.x renders import it); figure scripts move to this module one by one.
# Each script records STYLE_SHA256 in its numbers / style json.
import hashlib, io, re
import matplotlib as mpl
import matplotlib.transforms as mtrans
from matplotlib.font_manager import FontProperties, findfont

with open(__file__, "rb") as _fh:
    STYLE_SHA256 = hashlib.sha256(_fh.read()).hexdigest()
SOURCE_P4_SHA256 = "7a2e81a3db1537ed8bc9a9cfb608452e1024727d81eee816fb2abe174ac4146a"

MM = 1 / 25.4
WIDTH_FULL_MM, WIDTH_HALF_MM = 139.0, 67.0

# --- colours: Okabe-Ito; role colours fixed by the paper
BLACK, WHITE = "#000000", "#FFFFFF"
S3C, SPC, ZPC, VERM, SKY, PINK = "#E69F00", "#0072B2", "#009E73", "#D55E00", "#56B4E9", "#CC79A7"
MUTED = ["#7E8A94", "#9EA9B2", "#B7C0C7", "#8B949B", "#A6B0B7", "#C2CAD0", "#6F7A82", "#AEB7BD"]
MUTED_MARKERS = ["^", "v", "D", "P", "X", "<", ">", "h"]
GRID_ALPHA = 0.16
BAND_ALPHA = 0.14          # light background band for a regime (certified region, failure region)
BOX_FILL = "#F2F2F2"       # light grey fill for grouped regions in the block diagram
ROLE = {"S3": dict(color=S3C, marker="o"), "SPLIT6": dict(color=SPC, marker="s"), "Z+": dict(color=ZPC, marker="^")}

# --- type sizes (pt)
FS_TEXT, FS_TICK, FS_LEGEND = 8, 7.5, 7
FS_MIN = 7

# --- emphasis (P4): the series the argument is about is drawn at HERO weight, baselines recede
HERO_KW = dict(linewidth=1.7, markersize=4.6, markeredgewidth=0.5, markeredgecolor=WHITE, zorder=6)
MUTED_KW = dict(linewidth=0.75, markersize=2.6, alpha=0.9, zorder=2)
LW_AXIS, LW_REF = 0.6, 0.75
DASH_REF = (0, (3, 2))


def apply():
    """Global rcParams (P4 ieee_style.use_ieee with the BJ changes above).  Raises if the serif font does not resolve to Times New Roman."""
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "Nimbus Roman No9 L", "Liberation Serif"], "mathtext.fontset": "stix",
        "font.size": FS_TEXT, "axes.labelsize": FS_TEXT, "axes.titlesize": FS_TEXT, "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND, "figure.titlesize": FS_TEXT,
        "text.color": BLACK, "axes.labelcolor": BLACK, "xtick.color": BLACK, "ytick.color": BLACK, "axes.edgecolor": BLACK,
        "axes.linewidth": LW_AXIS, "lines.linewidth": 1.0, "lines.markersize": 3.2, "lines.markeredgewidth": 0.4, "patch.linewidth": 0.6,
        "xtick.major.width": LW_AXIS, "ytick.major.width": LW_AXIS, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
        "xtick.direction": "in", "ytick.direction": "in",
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False, "grid.linewidth": 0.4, "grid.alpha": GRID_ALPHA,
        "axes.titleweight": "normal", "axes.titlepad": 3,
        "legend.frameon": False, "legend.handlelength": 1.6, "legend.handletextpad": 0.5, "legend.columnspacing": 1.0, "legend.borderaxespad": 0.3,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600, "figure.dpi": 150,
        "savefig.bbox": None, "savefig.pad_inches": 0.0, "figure.facecolor": WHITE, "savefig.facecolor": WHITE,
        "figure.constrained_layout.use": True, "figure.constrained_layout.h_pad": 0.02, "figure.constrained_layout.w_pad": 0.02,
        "figure.constrained_layout.hspace": 0.03, "figure.constrained_layout.wspace": 0.03})
    f = findfont(FontProperties(family="serif"), fallback_to_default=False)
    if "Times New Roman" not in f: raise RuntimeError(f"serif font resolves to {f}, not Times New Roman")
    return f


def figsize(height_mm, width="full"):
    """BJ widths only: 'full' = 139 mm, 'half' = 67 mm."""
    if width not in ("full", "half"): raise ValueError(f"width must be 'full' or 'half', got {width!r}")
    return ((WIDTH_FULL_MM if width == "full" else WIDTH_HALF_MM) * MM, height_mm * MM)


def panel_label(ax, letter, title=None, x=0.0, y=1.0, pad_pt=3.0):
    """IEEE panel label '(a)' at the top left of the axes, optional short title after it on the same baseline, both 8 pt regular."""
    t = ax.annotate(f"({letter})", xy=(x, y), xycoords=ax.transAxes, xytext=(0, pad_pt), textcoords="offset points", fontsize=FS_TEXT, color=BLACK,
                    ha="left", va="bottom", annotation_clip=False)
    if title:
        ax.annotate(title, xy=(1, 0), xycoords=t, xytext=(3, 0), textcoords="offset points", fontsize=FS_TEXT, color=BLACK, ha="left", va="bottom", annotation_clip=False)
    return t


def hero(color, marker):
    """Kwargs for a HERO series (S3, SPLIT6, Z+): saturated colour, lw 1.7, marker with a white edge, top z-order."""
    return dict(color=color, marker=marker, **HERO_KW)


def muted(i, marker=None):
    """Kwargs for the i-th baseline / control: grey, lw 0.75, its own marker."""
    return dict(color=MUTED[i % len(MUTED)], marker=marker or MUTED_MARKERS[i % len(MUTED_MARKERS)], **MUTED_KW)


def band(ax, lo, hi, axis="y", color=SKY, alpha=BAND_ALPHA, **kw):
    """Light background band for a regime (e.g. the certified region above 1 - delta, or the failure region below it)."""
    f = ax.axhspan if axis == "y" else ax.axvspan
    return f(lo, hi, color=color, alpha=alpha, lw=0, zorder=0, **kw)


def ref_line(ax, value, axis="y", color=BLACK, dashed=True, lw=LW_REF):
    f = ax.axhline if axis == "y" else ax.axvline
    return f(value, color=color, lw=lw, ls=DASH_REF if dashed else "-", zorder=1)


def ygrid(ax):
    """The only grid allowed: y direction, alpha 0.16, behind the data."""
    ax.yaxis.grid(True, alpha=GRID_ALPHA, lw=0.4, color=BLACK); ax.set_axisbelow(True)


def mediabox_mm(pdf_bytes):
    """Width and height of the first /MediaBox in a PDF, in mm."""
    m = re.search(rb"/MediaBox\s*\[\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\]", pdf_bytes)
    if not m: raise RuntimeError("no /MediaBox in PDF")
    x0, y0, x1, y1 = (float(v) for v in m.groups())
    return (x1 - x0) / 72 * 25.4, (y1 - y0) / 72 * 25.4


def render(fig, width="full"):
    """PDF and 600 dpi PNG bytes; raises if the PDF /MediaBox width differs from the BJ width by more than 0.05 mm.  The caller writes them in 'x' mode."""
    if width not in ("full", "half"): raise ValueError(f"width must be 'full' or 'half', got {width!r}")
    buf = {}
    for ext, kw in (("pdf", {}), ("png", {"dpi": 600})):
        b = io.BytesIO(); fig.savefig(b, format=ext, **kw); buf[ext] = b.getvalue()
    w, h = mediabox_mm(buf["pdf"]); want = WIDTH_FULL_MM if width == "full" else WIDTH_HALF_MM
    if abs(w - want) > 0.05: raise RuntimeError(f"PDF width {w:.2f} mm, expected {want} mm")
    buf["mediabox_mm"] = (round(w, 2), round(h, 2))
    return buf
