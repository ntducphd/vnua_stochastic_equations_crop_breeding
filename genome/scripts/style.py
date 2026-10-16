"""Shared figure style (Nature-family specification for print; a slide preset for presentations).

Two output targets share the colours below:
  * print (setup / save_png): Figs 1-9 for journals, described next
  * slide (setup_slide / save_slide): 12 x 6 in canvas for 16:9 slides, posters and reports (step 14)

Print style, one for every figure:
  * canvas 183 mm (7.2 in) wide = Nature double column; height set per figure (<= 220 mm)
  * Arial (Helvetica / DejaVu Sans fallback), also for math text, so italic gene symbols and
    F_ST / H_e are set in Arial; text 5-7 pt at final size, panel labels 8 pt bold lowercase
  * axis lines 0.5 pt, top/right spines removed, small outward ticks, no grid
  * PNG, RGB, 600 dpi, white background, fixed canvas (exact 183 mm width)

Colour meaning (kept identical across figures):
  * CATEGORIES - validated categorical slots 1-5 (light mode; adjacent-pair CVD dE >= 9.1). Slots
    3-5 sit below 3:1 contrast, so every coloured mark also carries a category-specific marker
    shape and a black text label. Trait categories in Figs 1, 2 and 5; Fig 6 reuses the same slots
    for its five evidence tiers (defined in its own key).
  * DATA_GREY / DATA_GREY_LIGHT - neutral marks that carry no category (bars, dots, scan points).
  * DIV_NEG / DIV_MID / DIV_POS - diverging blue - neutral grey - red, centred on zero.
  * STRAND_SAME / STRAND_OPP - dark grey vs orange (also separable in greyscale).
Text is always black or grey, never coloured by series.
"""
import io
import re
from pathlib import Path

import matplotlib as mpl
from matplotlib.transforms import ScaledTranslation

CATEGORIES = {
    "1_starch_sugar": ("Starch & sugar (sweet / waxy)", "#2a78d6", "o"),
    "2_nutrition": ("Nutritional quality", "#eb6834", "s"),
    "3_pigment": ("Anthocyanin / pericarp colour", "#1baf7a", "D"),
    "4_plant_domestication": ("Plant type, flowering, domestication", "#eda100", "^"),
    "5_disease": ("Disease resistance", "#e87ba4", "v"),
}
INK = "#0b0b0b"          # primary text
INK2 = "#52514e"         # secondary text, axis lines
MUTED = "#8a8985"        # leader lines, reference lines
CHR_FILL = "#ecebe7"     # chromosome bodies
CHR_EDGE = "#a3a29c"
SURFACE = "#ffffff"
DATA_GREY = "#6b6a66"
DATA_GREY_LIGHT = "#b4b3ad"
DIV_NEG, DIV_MID, DIV_POS = "#2a78d6", "#f0efec", "#e34948"
STRAND_SAME, STRAND_OPP = "#52514e", "#eb6834"

# final-size type scale (pt)
FS = 6          # axis titles, legends, labels
FS_TICK = 5.5   # tick labels, in-plot annotations
FS_SMALL = 5    # dense annotations (gene labels in small multiples, coverage values)
FS_HEAD = 7     # chromosome / column headers
FS_PANEL = 8    # panel letters (bold, lowercase)
LW_AXIS = 0.5
WIDTH_IN = 7.2  # 183 mm
DPI = 600
MM = 25.4


def setup():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        # mathtext sets digits inside \mathit upright; the sf slot is mapped to Arial Italic so
        # italic_gene() can set whole symbols (sh2, ZmRap2.7, sweet4c) in italics
        "mathtext.sf": "Arial:italic",
        "font.size": FS,
        "text.color": INK,
        "axes.labelsize": FS,
        "axes.labelcolor": INK,
        "axes.labelpad": 2,
        "axes.titlesize": FS_HEAD,
        "axes.titlepad": 3,
        "axes.edgecolor": INK2,
        "axes.linewidth": LW_AXIS,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.facecolor": SURFACE,
        "axes.unicode_minus": True,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 2,
        "ytick.major.size": 2,
        "xtick.minor.size": 1.2,
        "ytick.minor.size": 1.2,
        "xtick.major.width": LW_AXIS,
        "ytick.major.width": LW_AXIS,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.pad": 1.5,
        "ytick.major.pad": 1.5,
        "lines.linewidth": 0.75,
        "patch.linewidth": LW_AXIS,
        "legend.fontsize": FS,
        "legend.title_fontsize": FS,
        "legend.frameon": False,
        "legend.handlelength": 1.4,
        "legend.handletextpad": 0.4,
        "legend.labelspacing": 0.35,
        "legend.columnspacing": 1.2,
        "legend.borderaxespad": 0.2,
        "legend.borderpad": 0.2,
        "image.interpolation": "nearest",
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "savefig.dpi": DPI,
    })


def italic_gene(label):
    """Gene symbols in italics (maize nomenclature), digits included; parentheses stay upright.
    Uses the math sf slot, which setup() maps to Arial Italic."""
    def it(tok):
        return r"$\mathsf{" + tok.replace(".", r".") + "}$"
    parts = re.split(r"(\s|\(|\))", label)
    return "".join(it(p) if p.strip() and p not in "()" else p for p in parts)


def chrom_number(name):
    m = re.fullmatch(r"(?i)chr0*(\d+)", str(name))
    return int(m.group(1)) if m else None


def minus(s):
    """Typographic minus sign for numbers written as text."""
    return str(s).replace("-", "−")


def spread(ys, gap, lo=float("-inf"), hi=float("inf")):
    """Push label positions apart by at least `gap`, keeping order and (if possible) bounds."""
    ys = list(ys)
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + gap)
    if ys and ys[-1] > hi:
        ys[-1] = hi
        for i in range(len(ys) - 2, -1, -1):
            ys[i] = min(ys[i], ys[i + 1] - gap)
    if ys and ys[0] < lo:
        ys[0] = lo
        for i in range(1, len(ys)):
            ys[i] = max(ys[i], ys[i - 1] + gap)
    return ys


def axes_in(fig, x, y, w, h, **kw):
    """Add axes positioned in inches from the bottom-left corner of the figure."""
    W, H = fig.get_size_inches()
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def axes_size_in(ax):
    """(width, height) of an axes in inches."""
    bb = ax.get_position()
    W, H = ax.figure.get_size_inches()
    return bb.width * W, bb.height * H


def panel_label(ax, letter, dx=-0.3, dy=0.0):
    """Bold 8-pt lowercase panel letter; its top-left corner sits dx, dy inches from the
    top-left corner of `ax` (negative dx = further left)."""
    tr = ax.transAxes + ScaledTranslation(dx, dy, ax.figure.dpi_scale_trans)
    ax.text(0, 1, letter, transform=tr, fontsize=FS_PANEL, fontweight="bold", color=INK,
            ha="left", va="top", clip_on=False)


def despine(ax, keep=("left", "bottom")):
    for side in ("left", "right", "top", "bottom"):
        ax.spines[side].set_visible(side in keep)


def save_png(fig, path):
    """Write `path` (.png) as RGB, 600 dpi, white background, on the fixed figure canvas.
    Prints the blank margin on each side (from the inked pixels) and flags artists that fall
    outside the canvas."""
    from PIL import Image, ImageChops
    path = Path(path).with_suffix(".png")
    fig.canvas.draw()
    W, H = fig.get_size_inches()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    clipped = bb.x0 < -0.005 or bb.y0 < -0.005 or bb.x1 > W + 0.005 or bb.y1 > H + 0.005
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=DPI, facecolor=SURFACE)
    buf.seek(0)
    im = Image.open(buf)
    im.load()
    if im.mode != "RGB":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.getchannel("A") if "A" in im.getbands() else None)
        im = bg
    im.save(path, format="PNG", dpi=(DPI, DPI), optimize=True)
    x0, y0, x1, y1 = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).getbbox()
    px = MM / DPI
    print(f"wrote {path.name}: {im.size[0]}x{im.size[1]} px = {W * MM:.0f} x {H * MM:.0f} mm; blank margins mm "
          f"left {x0 * px:.1f}, right {(im.size[0] - x1) * px:.1f}, top {y0 * px:.1f}, bottom {(im.size[1] - y1) * px:.1f}"
          + ("  <-- ARTIST OUTSIDE CANVAS" if clipped else ""))
    return path


# ---------------------------------------------------------------- presentation preset (step 14)
# A 12 x 6 in canvas fills a 13.33 x 7.5 in (16:9) slide under a title at 100 % scale. Type sizes are
# for that size: nothing below 12 pt, ticks 14 pt, axis titles 16 pt, headers 18 pt (readable from the
# back of a lecture room); lines >= 1 pt; markers >= 6 pt. Same colours and marker shapes as print.
SLIDE_W, SLIDE_H = 12.0, 6.0
SLIDE_DPI = 300
SL_FS, SL_TICK, SL_SMALL, SL_HEAD = 16, 14, 12, 18
SL_LW = 1.0
SL_MIN_PT = 12


def setup_slide():
    setup()
    mpl.rcParams.update({
        "font.size": SL_FS, "axes.labelsize": SL_FS, "axes.titlesize": SL_HEAD, "axes.labelpad": 5,
        "axes.titlepad": 6, "axes.linewidth": SL_LW, "xtick.labelsize": SL_TICK, "ytick.labelsize": SL_TICK,
        "xtick.major.size": 5, "ytick.major.size": 5, "xtick.minor.size": 3, "ytick.minor.size": 3,
        "xtick.major.width": SL_LW, "ytick.major.width": SL_LW, "xtick.minor.width": 0.8,
        "ytick.minor.width": 0.8, "xtick.major.pad": 4, "ytick.major.pad": 4, "lines.linewidth": 1.8,
        "patch.linewidth": SL_LW, "legend.fontsize": SL_TICK, "legend.title_fontsize": SL_TICK,
        "legend.handlelength": 1.6, "legend.handletextpad": 0.5, "legend.labelspacing": 0.5,
        "savefig.dpi": SLIDE_DPI,
    })


def save_slide(fig, path, min_pt=SL_MIN_PT):
    """Write `path` (.png) at 300 dpi on the fixed 12 x 6 in canvas. Reports the smallest font size used and
    flags any text below `min_pt` or any artist outside the canvas."""
    from PIL import Image
    path = Path(path).with_suffix(".png")
    fig.canvas.draw()
    W, H = fig.get_size_inches()
    r = fig.canvas.get_renderer()
    bb = fig.get_tightbbox(r)
    clipped = bb.x0 < -0.005 or bb.y0 < -0.005 or bb.x1 > W + 0.005 or bb.y1 > H + 0.005
    sizes = [t.get_fontsize() for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text().strip()]
    small = min(sizes) if sizes else None
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=SLIDE_DPI, facecolor=SURFACE)
    buf.seek(0)
    im = Image.open(buf)
    im.load()
    im.convert("RGB").save(path, format="PNG", dpi=(SLIDE_DPI, SLIDE_DPI), optimize=True)
    print(f"wrote {path.name}: {im.size[0]}x{im.size[1]} px = {W:.1f} x {H:.1f} in; smallest text {small:.1f} pt"
          + ("  <-- TEXT BELOW MINIMUM" if small is not None and small < min_pt - 0.01 else "")
          + ("  <-- ARTIST OUTSIDE CANVAS" if clipped else ""))
    return path
