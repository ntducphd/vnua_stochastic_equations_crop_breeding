# -*- coding: utf-8 -*-
"""Model figures of the review.

  R01_fig_drift_selection  a, stationary density of the Wright-Fisher diffusion with mutation for
                           4 Ne s = 0.5, 2 and 8 (Ne = 100, mu1 = mu2 = 0.001, so 4 Ne mu = 0.4);
                           b, variance of the allele frequency under that density against 4 Ne s;
                           c, Kimura's fixation probability from x0 = 0.5 against 4 Ne s
  R01_fig_wf_example       the Wright-Fisher worked example (Ne = 100, s = 0.15, x0 = 0.5):
                           a, 30 Euler-Maruyama paths of the SDE (dt = 0.1 generation, seed 42) and their mean;
                           b, exact distribution of the discrete chain (2 Ne = 200 gene copies, genic
                           selection) at generations 2, 4, 6 and 10; c, Kimura's fixation probability
                           over x0 and Ne at s = 0.03
  R01_fig_family_map       schematic of the five families of stochastic equations, what each answers
                           in breeding, and how they are related

Every figure is drawn with fig_style (160 mm wide, Arial 7 pt, panel letters 8 pt bold, Okabe-Ito
colours) and saved to ../figures/ as a vector PDF and a 600 dpi PNG at the drawn size. The script
prints the values that the captions quote and writes the package versions to ../figures/ENVIRONMENT.txt.
The only random numbers are the Euler-Maruyama increments of R01_fig_wf_example a (seed 42).

Usage: python -I -B scripts/generate_figures.py
"""
import os
import platform
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402
from scipy.integrate import quad  # noqa: E402
from scipy.special import betaln, hyp1f1  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # python -I does not put the script folder on the path
from fig_style import GREY, INK, OI, PT, apply, new_figure, panel, save  # noqa: E402
from wf_exact import transition as wf_transition  # noqa: E402  (exact discrete chain)

OUTDIR = os.path.normpath(os.path.join(HERE, "..", "figures"))
SEED_PATHS = 42


# ---------------------------------------------------------------------------------------------
# Wright-Fisher diffusion: stationary density, its moments, Kimura's fixation probability
# ---------------------------------------------------------------------------------------------
NE_STAT = 100      # R01_fig_drift_selection a, b: Ne fixed and s varied, so 4 Ne mu = 0.4 in every curve
MU_STAT = 0.001


def wf_stationary(p_arr, Nes4, Ne=NE_STAT, mu1=MU_STAT, mu2=MU_STAT):
    """Normalised stationary density
    phi(x) = x^(a - 1) (1 - x)^(b - 1) exp(4 Ne s x) / [B(a, b) 1F1(a; a + b; 4 Ne s)],
    a = 4 Ne mu1, b = 4 Ne mu2 (mu1: mutation towards the allele whose frequency is x)."""
    a, b = 4 * Ne * mu1, 4 * Ne * mu2
    log_norm = betaln(a, b) + np.log(hyp1f1(a, a + b, Nes4))
    return np.exp((a - 1) * np.log(p_arr) + (b - 1) * np.log1p(-p_arr) + Nes4 * p_arr - log_norm)


def wf_moments(Nes4, Ne=NE_STAT, mu1=MU_STAT, mu2=MU_STAT):
    """Mean and variance of the stationary density, from ratios of confluent hypergeometric functions."""
    a, b = 4 * Ne * mu1, 4 * Ne * mu2
    f0 = hyp1f1(a, a + b, Nes4)
    m1 = a / (a + b) * hyp1f1(a + 1, a + b + 1, Nes4) / f0
    m2 = a * (a + 1) / ((a + b) * (a + b + 1)) * hyp1f1(a + 2, a + b + 2, Nes4) / f0
    return m1, m2 - m1 ** 2


def wf_tail_masses(Nes4, Ne=NE_STAT, mu1=MU_STAT, mu2=MU_STAT, edge=0.1):
    """P(x < edge) and P(x > 1 - edge) under the stationary density. The endpoint singularities
    x^(a - 1) and (1 - x)^(b - 1) are passed to quad as algebraic weights."""
    a, b = 4 * Ne * mu1, 4 * Ne * mu2
    log_norm = betaln(a, b) + np.log(hyp1f1(a, a + b, Nes4))
    low = quad(lambda x: np.exp((b - 1) * np.log1p(-x) + Nes4 * x - log_norm), 0, edge,
               weight="alg", wvar=(a - 1, 0))[0]
    high = quad(lambda x: np.exp((a - 1) * np.log(x) + Nes4 * x - log_norm), 1 - edge, 1,
                weight="alg", wvar=(0, b - 1))[0]
    return low, high


def kimura_pfix(x0, Nes4_arr):
    """Kimura's fixation probability (1 - exp(-4 Ne s x0)) / (1 - exp(-4 Ne s)) for each 4 Ne s."""
    result = np.zeros_like(Nes4_arr, dtype=float)
    for i, g in enumerate(Nes4_arr):
        if abs(g) < 1e-8:
            result[i] = x0
        else:
            num = 1 - np.exp(-g * x0)
            den = 1 - np.exp(-g)
            result[i] = num / den if abs(den) > 1e-12 else x0
    return result


def euler_maruyama_wf(s, Ne, x0, n_gen, n_rep, dt=1.0, seed=SEED_PATHS):
    """Euler-Maruyama paths of dx = s x (1 - x) dt + sqrt(x (1 - x) / (2 Ne)) dW, kept in [0, 1]."""
    rng = np.random.default_rng(seed)
    n_steps = int(n_gen / dt)
    paths = np.zeros((n_rep, n_steps + 1))
    paths[:, 0] = x0
    for t in range(n_steps):
        x = np.clip(paths[:, t], 1e-6, 1 - 1e-6)
        drift = s * x * (1 - x) * dt
        diff = np.sqrt(x * (1 - x) / (2 * Ne)) * np.sqrt(dt)
        paths[:, t + 1] = np.clip(x + drift + diff * rng.normal(0, 1, n_rep), 0, 1)
    return paths


def _letters(fig, axes, letters):
    """Panel letters above the top left of each panel, flush with the panel's leftmost text
    (axis label or tick labels), so that they never touch a tick label."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for ax, letter in zip(axes, letters):
        probe = ax.text(0, 0, letter, fontsize=8, fontweight="bold")
        w_px = probe.get_window_extent(r).width
        probe.remove()
        x_left = ax.transAxes.inverted().transform((ax.get_tightbbox(r).x0, 0))[0]
        width_px = ax.get_window_extent(r).width
        panel(ax, letter, x=x_left + w_px / width_px, y=1.0)


# ---------------------------------------------------------------------------------------------
# R01_fig_drift_selection
# ---------------------------------------------------------------------------------------------
def fig_drift_selection():
    fig, axes = new_figure(nrows=1, ncols=3, height_mm=55)
    p_inner = np.linspace(0.0005, 0.9995, 2000)

    # a: stationary densities; colour and line style both distinguish the curves
    ax = axes[0]
    regimes = [(0.5, OI["blue"], "solid", r"$4N_e s$ = 0.5"),
               (2.0, OI["orange"], "dashed", r"$4N_e s$ = 2"),
               (8.0, OI["vermilion"], "dotted", r"$4N_e s$ = 8")]
    pi_max, pi_min = 0, np.inf
    for Nes4, col, ls, lbl in regimes:
        dens = wf_stationary(p_inner, Nes4)
        ax.plot(p_inner, dens, color=col, ls=ls, label=lbl)
        pi_max, pi_min = max(pi_max, dens.max()), min(pi_min, dens.min())
    ax.set_yscale("log")
    ax.set_ylim(pi_min / 2, pi_max * 5)
    ax.set_xlim(0, 1)
    ax.set_xlabel(r"Allele frequency $x$")
    ax.set_ylabel("Stationary density")
    ax.legend(loc="upper center")

    # b: variance of x under the stationary density (analytic moments)
    ax = axes[1]
    Nes4_range = np.logspace(-1, 1, 400)
    variances = np.array([wf_moments(g)[1] for g in Nes4_range])
    slope = np.gradient(variances, np.log(Nes4_range))
    steep = Nes4_range[np.argmin(slope)]
    ax.plot(Nes4_range, variances, color=INK)
    ax.axvline(x=1.0, color=GREY, lw=0.6, ls="--")
    ax.text(0.88, variances.max() * 0.30, "conventional\nboundary", ha="right", fontstyle="italic")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    ax.set_ylim(0, variances.max() * 1.15)
    ax.set_xlabel(r"$4N_e s$")
    ax.set_ylabel(r"Var($x$)")

    # c: Kimura's fixation probability from x0 = 0.5
    ax = axes[2]
    Nes4_c = np.linspace(-4, 8, 300)
    ax.plot(Nes4_c, kimura_pfix(0.5, Nes4_c), color=OI["blue"])
    ax.axhline(0.5, color=GREY, lw=0.6, ls=":")
    ax.axvline(0, color=GREY, lw=0.6, ls=":")
    ax.set_xticks([-4, 0, 4, 8])
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$4N_e s$")
    ax.set_ylabel("Fixation probability")

    _letters(fig, axes, "abc")
    path = save(fig, OUTDIR, "R01_fig_drift_selection")
    print("R01_fig_drift_selection (%s)" % path)
    print("  a: Ne = %d, mu1 = mu2 = %g, 4 Ne mu = %.1f" % (NE_STAT, MU_STAT, 4 * NE_STAT * MU_STAT))
    for Nes4, _, _, _ in regimes:
        low, high = wf_tail_masses(Nes4)
        m, v = wf_moments(Nes4)
        print("     4Nes = %-3g P(x < 0.1) = %.3f, P(x > 0.9) = %.3f, mean %.3f, variance %.3f"
              % (Nes4, low, high, m, v))
    print("  b: the variance falls most steeply at 4Nes = %.2f; Var(x) at 4Nes = 0.1, 1, 2, 3, 5: %s"
          % (steep, ", ".join("%.3f" % wf_moments(g)[1] for g in (0.1, 1, 2, 3, 5))))
    print("  c: P_fix from x0 = 0.5 at 4Nes = 0.5, 1, 2, 8: %s"
          % ", ".join("%.3f" % v for v in kimura_pfix(0.5, np.array([0.5, 1.0, 2.0, 8.0]))))


# ---------------------------------------------------------------------------------------------
# R01_fig_wf_example
# ---------------------------------------------------------------------------------------------
def fig_wf_example():
    fig, axes = new_figure(nrows=1, ncols=3, height_mm=58)
    s, Ne, x0, n_gen, n_rep = 0.15, 100, 0.5, 10, 30

    # a: Euler-Maruyama paths and their mean
    ax = axes[0]
    paths = euler_maruyama_wf(s, Ne, x0, n_gen, n_rep, dt=0.1)
    t_arr = np.linspace(0, n_gen, paths.shape[1])
    for i in range(n_rep):
        ax.plot(t_arr, paths[i], color=GREY, lw=0.4, alpha=0.6)
    mean_path = paths.mean(axis=0)
    ax.plot(t_arr, mean_path, color=OI["vermilion"], lw=1.4, label="Mean path")
    ax.axhline(0.95, color=GREY, lw=0.6, ls="--")
    ax.text(0.3, 0.964, "target frequency", transform=ax.get_yaxis_transform(), ha="left", va="bottom")
    ax.set_xlim(0, n_gen)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Generation")
    ax.set_ylabel(r"Allele frequency $x$")
    ax.legend(loc="lower right")

    # b: exact distribution of the discrete chain (no random numbers)
    ax = axes[1]
    snapshots = [(2, OI["orange"]), (4, OI["green"]), (6, OI["blue"]), (10, OI["black"])]
    P = wf_transition(Ne, s)                       # 2 Ne + 1 = 201 states, genic selection
    two_n = 2 * Ne
    p_grid = np.arange(two_n + 1) / two_n
    v = np.zeros(two_n + 1)
    v[int(round(x0 * two_n))] = 1.0
    t_now = 0
    chain_mean = {}
    for t_snap, col in snapshots:
        while t_now < t_snap:
            v = v @ P
            t_now += 1
        chain_mean[t_snap] = float(v @ p_grid)
        density = v * two_n                        # probability per unit allele frequency
        ax.plot(p_grid, density, color=col, label=r"$t$ = %d" % t_snap)
        ax.fill_between(p_grid, density, alpha=0.10, color=col, lw=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"Allele frequency $x$")
    ax.set_ylabel("Probability density")
    ax.legend(loc="upper left")

    # c: Kimura's fixation probability over x0 and Ne (log axis) at s = 0.03
    ax = axes[2]
    x0_arr = np.linspace(0.05, 0.95, 60)
    Ne_fine = np.logspace(np.log10(10), np.log10(500), 40)
    s_c = 0.03
    pfix_mat = np.zeros((len(Ne_fine), len(x0_arr)))
    for i, Ne_i in enumerate(Ne_fine):
        Nes4 = 4 * Ne_i * s_c
        for j, x0_j in enumerate(x0_arr):
            if abs(Nes4) < 1e-8:
                pfix_mat[i, j] = x0_j
            else:
                num = 1 - np.exp(-Nes4 * x0_j)
                den = 1 - np.exp(-Nes4)
                pfix_mat[i, j] = num / den if abs(den) > 1e-12 else x0_j
    im = ax.pcolormesh(x0_arr, Ne_fine, pfix_mat, cmap="viridis", vmin=0, vmax=1, shading="auto",
                       rasterized=True)
    ax.set_yscale("log")
    ax.set_yticks([10, 20, 50, 100, 200, 500])
    ax.set_yticklabels(["10", "20", "50", "100", "200", "500"])
    half = np.diff(np.log10(Ne_fine))[0] / 2      # mesh cells reach half a step beyond the grid
    ax.set_ylim(Ne_fine[0] * 10 ** -half, Ne_fine[-1] * 10 ** half)
    ax.text(0.04, 0.04, r"$s$ = 0.03", transform=ax.transAxes, color="white", ha="left", va="bottom")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Fixation probability")
    cbar.outline.set_linewidth(0.6)
    ax.set_xlabel(r"Initial frequency $x_0$")
    ax.set_ylabel(r"Effective population size $N_e$")

    _letters(fig, axes, "abc")
    path = save(fig, OUTDIR, "R01_fig_wf_example")
    print("R01_fig_wf_example (%s)" % path)
    print("  a: Ne = %d, s = %.2f, x0 = %.1f, %d Euler-Maruyama paths, dt = 0.1 generation, seed %d; "
          "mean of the paths at generation %d = %.3f" % (Ne, s, x0, n_rep, SEED_PATHS, n_gen, mean_path[-1]))
    print("  b: exact chain mean at generations 2, 4, 6, 10: %s"
          % ", ".join("%.3f" % chain_mean[t] for t, _ in snapshots))
    k10, k50 = kimura_pfix(0.05, np.array([4 * 10 * s, 4 * 50 * s]))
    print("  c: drawn at s = %.2f; at the example's s = %.2f and x0 = 0.05, P_fix = %.3f (Ne = 10) and %.3f (Ne = 50)"
          % (s_c, s, k10, k50))


# ---------------------------------------------------------------------------------------------
# R01_fig_family_map
# ---------------------------------------------------------------------------------------------
# Each family: title, equations or models (regular), what it answers in breeding (italic).
# Colours follow the family colours of the evidence map (13_evidence_map.py, fig_style.CYCLE order):
# master equations blue, diffusion vermilion, branching orange, control purple, simulation sky blue.
FAMILIES = [
    ("Master equations and Markov chains", OI["blue"],
     ["Wright–Fisher and Moran chains; selfing and doubled-haploid",
      "chains; chemical master equation"],
     ["Exact finite-population probabilities"]),
    ("Diffusion approximations", OI["vermilion"],
     ["Wright–Fisher SDE; Fokker–Planck (forward) and",
      "Kolmogorov (backward) equations"],
     ["Fixation probabilities, absorption times, response and its spread"]),
    ("Branching processes", OI["orange"],
     ["Rare alleles, edited events and introgressed segments",
      "while they are rare"],
     [r"Establishment probability $\pi\;\approx\;2s$"]),
    ("Stochastic optimal control", OI["purple"],
     ["Bellman equation in discrete generations;",
      "HJB equation in the diffusion limit"],
     ["Multi-generation selection, mating and resource allocation"]),
]
SIMULATION = (["Stochastic simulation", "of programmes"], OI["sky"],
              ["Individual-based forward", "simulation; numerical", "counterpart of all four"],
              ["Whole-programme designs", "that have no closed form"])

LINE = 3.05      # mm between text lines at 7 pt
PAD = 2.2        # mm between a box edge and its first or last line


def _tint(colour, f=0.12):
    r, g, b = mcolors.to_rgb(colour)
    return (1 - f + f * r, 1 - f + f * g, 1 - f + f * b)


def _box(ax, x, y, w, h, colour, title, body, answers):
    """Rounded box with a bold title, regular body lines and italic answer lines, centred."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.5",
                                facecolor=_tint(colour), edgecolor=colour, linewidth=0.9))
    lines = [(t, "bold") for t in title] + [(t, "normal") for t in body] + [(t, "italic") for t in answers]
    n = len(lines)
    y_first = y + h / 2 + (n - 1) * LINE / 2
    texts = []
    for k, (txt, style) in enumerate(lines):
        kw = {"fontweight": "bold"} if style == "bold" else ({"fontstyle": "italic"} if style == "italic" else {})
        texts.append(ax.text(x + w / 2, y_first - k * LINE, txt, ha="center", va="center", **kw))
    return texts


def _arrow(ax, p0, p1, arrows, colour=INK):
    """Straight arrow from p0 to p1 (mm); its end points are added to `arrows` for check_layout."""
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=7, linewidth=0.8, color=colour,
                                 shrinkA=0, shrinkB=0, joinstyle="miter", capstyle="butt"))
    arrows.append((p0, p1))


def _segments(ax, arrows):
    """Straight segments of every line and arrow, in data coordinates (for the overlap check)."""
    segs = list(arrows)
    for ln in ax.lines:
        xs, ys = ln.get_data()
        segs += [((xs[i], ys[i]), (xs[i + 1], ys[i + 1])) for i in range(len(xs) - 1)]
    return segs


def _segment_hits_box(seg, box):
    """True if the segment passes through the rectangle (x0, y0, x1, y1); sampled densely."""
    (xa, ya), (xb, yb) = seg
    t = np.linspace(0, 1, 400)
    xs, ys = xa + t * (xb - xa), ya + t * (yb - ya)
    x0, y0, x1, y1 = box
    return bool(np.any((xs > x0) & (xs < x1) & (ys > y0) & (ys < y1)))


def check_layout(fig, ax, owned, arrows):
    """Report text that overlaps other text, leaves its box, crosses a line or arrow, or leaves the
    canvas. `owned` maps each text inside a box to that box (x, y, w, h); other texts are labels."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    tb = []
    for t in ax.texts:
        bb = t.get_window_extent(r)
        (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
        tb.append((t, (x0, y0, x1, y1)))
    problems = []
    for i in range(len(tb)):
        for j in range(i + 1, len(tb)):
            a, b = tb[i][1], tb[j][1]
            if a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]:
                problems.append("text overlap: %r / %r" % (tb[i][0].get_text(), tb[j][0].get_text()))
    W, H = ax.get_xlim()[1], ax.get_ylim()[1]
    boxes = sorted(set(owned.values()))
    segs = _segments(ax, arrows)
    for t, (x0, y0, x1, y1) in tb:
        if x0 < 0 or y0 < 0 or x1 > W or y1 > H:
            problems.append("outside the canvas: %r" % t.get_text())
        if t in owned:
            bx, by, bw, bh = owned[t]
            if x0 < bx + 0.8 or x1 > bx + bw - 0.8 or y0 < by + 0.5 or y1 > by + bh - 0.5:
                problems.append("text leaves its box: %r" % t.get_text())
        else:
            for bx, by, bw, bh in boxes:
                if x0 < bx + bw and bx < x1 and y0 < by + bh and by < y1:
                    problems.append("label overlaps a box: %r" % t.get_text())
            pad = (x0 - 0.4, y0 - 0.4, x1 + 0.4, y1 + 0.4)
            if any(_segment_hits_box(sg, pad) for sg in segs):
                problems.append("label touches a line or arrow: %r" % t.get_text())
    return problems


def fig_family_map():
    W, H = 160.0, 88.0
    fig, ax = new_figure(height_mm=H, constrained_layout=False)
    ax.set_position([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    # four families stacked in the order of the review; simulation as one tall box to their right
    x0, w = 8.0, 92.0                       # stack
    h, top, bottom = 16.0, H - 1.0, 1.0
    gap = (top - bottom - 4 * h) / 3
    owned, mids, frames = {}, [], []
    for k, (title, colour, body, answers) in enumerate(FAMILIES):
        y = top - h - k * (h + gap)
        for t in _box(ax, x0, y, w, h, colour, [title], body, answers):
            owned[t] = (x0, y, w, h)
        mids.append(y + h / 2)
        frames.append((y, y + h))
    (yA, yB, yC, yD) = mids
    (A, B, C, D) = frames                   # (bottom, top) of each box
    xm = x0 + w / 2
    arrows = []
    # arrow labels upright; italic is kept for what each family answers
    lab = {"ha": "left", "va": "center", "color": INK}

    # chains -> diffusion
    _arrow(ax, (xm, A[0]), (xm, B[1]), arrows)
    ax.text(xm + 2, (A[0] + B[1]) / 2, r"large $N_e$, weak selection", **lab)
    # diffusion -> branching and branching -> diffusion, side by side
    _arrow(ax, (xm - 3, B[0]), (xm - 3, C[1]), arrows)
    ax.text(xm - 5, (B[0] + C[1]) / 2, "allele still rare", **dict(lab, ha="right"))
    _arrow(ax, (xm + 3, C[1]), (xm + 3, B[0]), arrows)
    ax.text(xm + 5, (B[0] + C[1]) / 2, r"agree for small $s$: $u\;\approx\;\pi\;\approx\;2s$", **lab)

    # chains and diffusion -> control: a bus along the left edge
    xb = 4.0
    ax.plot([x0, xb, xb], [yA, yA, yD], color=INK, lw=0.8, solid_joinstyle="miter")
    ax.plot([x0, xb], [yB, yB], color=INK, lw=0.8)
    _arrow(ax, (xb, yD), (x0, yD), arrows)
    ax.text(xb + 2, (C[0] + D[1]) / 2, "transition law of a chain or a diffusion", **lab)

    # all four -> simulation: one short grey arrow from each family into a box spanning the stack
    xs = x0 + w + 9.0
    ws = W - 1.0 - xs
    for ym in mids:
        _arrow(ax, (x0 + w, ym), (xs, ym), arrows, colour=GREY)
    title, colour, body, answers = SIMULATION
    ys, hs = D[0], A[1] - D[0]
    for t in _box(ax, xs, ys, ws, hs, colour, title, body, answers):
        owned[t] = (xs, ys, ws, hs)

    problems = check_layout(fig, ax, owned, arrows)
    path = save(fig, OUTDIR, "R01_fig_family_map")
    print("R01_fig_family_map (%s): %s" % (path, "layout check passed (no overlaps)" if not problems
                                           else "LAYOUT PROBLEMS:\n    " + "\n    ".join(problems)))
    return problems


def write_environment():
    with open(os.path.join(OUTDIR, "ENVIRONMENT.txt"), "w", encoding="utf-8") as fh:
        fh.write("generate_figures.py: R01_fig_drift_selection, R01_fig_wf_example, R01_fig_family_map\n")
        fh.write("python %s\nnumpy %s\nscipy %s\nmatplotlib %s\n" % (
            platform.python_version(), np.__version__, scipy.__version__, matplotlib.__version__))
        fh.write("random seed: %d (Euler-Maruyama paths of R01_fig_wf_example a); every other panel is "
                 "deterministic\n" % SEED_PATHS)


if __name__ == "__main__":
    apply()
    assert PT == 7
    fig_drift_selection()
    fig_wf_example()
    problems = fig_family_map()
    write_environment()
    print("PDF and 600 dpi PNG written to %s" % OUTDIR)
    sys.exit(1 if problems else 0)
