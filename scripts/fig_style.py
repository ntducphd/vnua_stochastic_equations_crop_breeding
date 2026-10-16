"""One style for every figure of the review.

Figures are drawn at the width they are printed, the 160 mm text column of an A4 page with
2.5 cm margins, so that no type is scaled after drawing:

  * Arial, 7 pt for labels, ticks and legends; panel letters 8 pt bold lowercase
  * L-shaped axes, 0.6 pt lines, ticks outward
  * the Okabe-Ito palette, which stays distinct under the common forms of colour blindness
  * vector PDF with the fonts embedded as TrueType, and a 600 dpi PNG for the Word file

Usage:
    from fig_style import apply, new_figure, panel, save, OI
    apply()
    fig, axes = new_figure(nrows=1, ncols=3, height_mm=55)
    ...
    panel(axes[0], "a")
    save(fig, outdir, "R01_fig_name")
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt

MM = 1 / 25.4                      # inches per millimetre
WIDTH_MM = 160.0                   # the text column
INK = "#000000"
GREY = "#7F7F7F"
LIGHT = "#D9D9D9"

# Okabe-Ito
OI = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermilion": "#D55E00",
    "purple": "#CC79A7",
}
CYCLE = [OI["blue"], OI["vermilion"], OI["green"], OI["orange"], OI["purple"], OI["sky"], OI["black"]]

PT = 7          # text
PT_LETTER = 8   # panel letters
LW = 0.6        # axes, ticks
LW_DATA = 1.0   # data lines


def apply():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": PT,
        "axes.titlesize": PT,
        "axes.labelsize": PT,
        "xtick.labelsize": PT,
        "ytick.labelsize": PT,
        "legend.fontsize": PT,
        "legend.title_fontsize": PT,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "axes.unicode_minus": True,
        "axes.linewidth": LW,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": mpl.cycler(color=CYCLE),
        "axes.titlepad": 3,
        "axes.labelpad": 2,
        "xtick.major.width": LW,
        "ytick.major.width": LW,
        "xtick.minor.width": LW * 0.8,
        "ytick.minor.width": LW * 0.8,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "lines.linewidth": LW_DATA,
        "lines.markersize": 3,
        "patch.linewidth": LW,
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.borderaxespad": 0.3,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,          # TrueType, embedded and editable
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def new_figure(nrows=1, ncols=1, height_mm=60.0, width_mm=WIDTH_MM, **kw):
    kw.setdefault("constrained_layout", True)
    fig, axes = plt.subplots(nrows, ncols, figsize=(width_mm * MM, height_mm * MM), **kw)
    return fig, axes


def panel(ax, letter, x=-0.02, y=1.0):
    """Bold lowercase panel letter at the top left, outside the axes."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=PT_LETTER, fontweight="bold",
            ha="right", va="bottom")


def save(fig, outdir, stem):
    """PDF (vector) and PNG (600 dpi) at the drawn size; no tight bounding box, so the printed
    width stays the 160 mm the figure was drawn at."""
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, stem + ".pdf"))
    fig.savefig(os.path.join(outdir, stem + ".png"), dpi=600)
    plt.close(fig)
    return os.path.join(outdir, stem + ".png")
