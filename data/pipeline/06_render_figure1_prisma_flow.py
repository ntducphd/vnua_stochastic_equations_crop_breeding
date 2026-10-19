"""
06_render_figure1_prisma_flow.py -- Figure 1: PRISMA 2020 flow diagram of the literature search (identification,
deduplication and title/abstract screening; full texts were not assessed).

Every count is computed from the data: identification from the Scopus and Web of Science exports in raw/, and
deduplication and screening from pipeline_output/03_final_decisions.csv, the decisions of the screening rules
(pipeline/02_screen_records.py; manuscript Section 2.2).

Drawn with ../scripts/fig_style.py at the printed width (160 mm, Arial 7 pt) and saved as vector PDF and 600 dpi PNG.

Usage: python -I pipeline/06_render_figure1_prisma_flow.py   (run from data/, after step 02)
"""
import glob
import os
import sys

import pandas as pd

sys.dont_write_bytecode = True   # no __pycache__ folder in the deposited scripts folder
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, "scripts"))
from fig_style import apply, new_figure, save, INK, LW  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

OUT_DIR = "../figures"
STEM = "R01_fig_prisma"

REASONS = {   # final_exclusion_reason -> wording in the figure, in display order
    "not_genetics_domain": "outside the genetics and breeding domain",
    "keyword_coincidence": "search terms matched incidentally",
    "not_stochastic_framework": "no stochastic model",
    "wrong_document_type": "not a research or review article",
}

# ---- counts ----
n_scopus = len(pd.read_csv(glob.glob("raw/scopus/*.csv")[0], encoding="utf-8-sig"))
n_wos = len(pd.read_excel(glob.glob("raw/WoS/*.xls")[0]))
final = pd.read_csv("pipeline_output/03_final_decisions.csv", encoding="utf-8-sig")
n_screened = len(final)
n_duplicates = n_scopus + n_wos - n_screened
decisions = final["final_decision"].value_counts()
n_retained, n_uncertain, n_excluded = (int(decisions.get(k, 0)) for k in ("Include", "Maybe", "Exclude"))
assert n_retained + n_uncertain + n_excluded == n_screened
reasons = final.loc[final["final_decision"] == "Exclude", "final_exclusion_reason"].value_counts()
assert set(reasons.index) <= set(REASONS), set(reasons.index) - set(REASONS)


def n(x):
    return f"(n = {x:,})"


# ---- layout (millimetres, y downwards) ----
apply()
W, H = 160.0, 80.0
fig, ax = new_figure(height_mm=H, constrained_layout=False)
ax.set_position([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.axis("off")

STAGE_X, STAGE_W = 2.0, 6.0          # stage bars at the left edge
MAIN_X, MAIN_W = 13.0, 66.0          # main column
SIDE_X, SIDE_W = 91.0, 67.0          # side column
MAIN_C, SIDE_C = MAIN_X + MAIN_W / 2, SIDE_X + SIDE_W / 2
STAGE_FILL = "#DDEFFA"               # pale tint of Okabe-Ito sky blue


def box(x, y, w, h, lines, dashed=False):
    ax.add_patch(Rectangle((x, y), w, h, facecolor="white", edgecolor=INK, linewidth=LW,
                           linestyle=(0, (3, 2)) if dashed else "solid", zorder=2))
    ax.text(x + w / 2, y + h / 2, "\n".join(lines), ha="center", va="center", linespacing=1.35, zorder=3)


def line(xs, ys, head=True):
    """Connector through the points (xs, ys); the last segment ends in an arrow head."""
    if len(xs) > 2:
        ax.plot(xs[:-1], ys[:-1], color=INK, lw=LW, solid_capstyle="butt", zorder=1)
    ax.annotate("", xy=(xs[-1], ys[-1]), xytext=(xs[-2], ys[-2]), zorder=1,
                arrowprops=dict(arrowstyle="-|>" if head else "-", color=INK, lw=LW, mutation_scale=7,
                                shrinkA=0, shrinkB=0))


def stage(y0, y1, label):
    ax.add_patch(FancyBboxPatch((STAGE_X, y0), STAGE_W, y1 - y0, boxstyle="round,pad=0,rounding_size=1.2",
                                facecolor=STAGE_FILL, edgecolor="none", zorder=1))
    ax.text(STAGE_X + STAGE_W / 2, (y0 + y1) / 2, label, rotation=90, ha="center", va="center", fontweight="bold")


# Identification: y 4-21
stage(4.0, 21.0, "Identification")
box(MAIN_X, 4.0, MAIN_W, 17.0, ["Records identified from:", f"Scopus {n(n_scopus)}", f"Web of Science {n(n_wos)}"])
box(SIDE_X, 6.0, SIDE_W, 13.0, ["Records removed before screening:", f"duplicate records {n(n_duplicates)}"])
line([MAIN_X + MAIN_W, SIDE_X], [12.5, 12.5])

# Screening: y 24.5-67
stage(24.5, 67.0, "Screening")
box(MAIN_X, 29.0, MAIN_W, 12.0, ["Records screened", n(n_screened)])
line([MAIN_C, MAIN_C], [21.0, 29.0])
box(SIDE_X, 24.5, SIDE_W, 21.0, [f"Records excluded {n(n_excluded)}:"]
    + [f"{REASONS[k]} {n(int(reasons.get(k, 0)))}" for k in REASONS])
line([MAIN_X + MAIN_W, SIDE_X], [35.0, 35.0])
box(MAIN_X, 55.0, MAIN_W, 12.0, ["Records retained as within scope", n(n_retained)])
line([MAIN_C, MAIN_C], [41.0, 55.0])
box(SIDE_X, 55.0, SIDE_W, 12.0, ["Uncertain", n(n_uncertain)], dashed=True)
line([MAIN_C, SIDE_C, SIDE_C], [49.5, 49.5, 55.0])

ax.text(MAIN_X, 73.5,
        "Screening used titles and abstracts. Full texts were not assessed, so the uncertain records were neither "
        "included nor excluded.", ha="left", va="center")

out = save(fig, OUT_DIR, STEM)
print(f"identified {n_scopus + n_wos:,} (Scopus {n_scopus:,}, Web of Science {n_wos:,}); duplicates removed "
      f"{n_duplicates:,}; screened {n_screened:,}; retained {n_retained:,}; uncertain {n_uncertain:,}; "
      f"excluded {n_excluded:,} {reasons.to_dict()}")
print(f"Saved {out} and the PDF")
