"""The a1-yz1-x1-sh2 interval on chromosome 3 in the 17 assembled genomes (purple supersweet linkage).

Why: the dominant A1 allele (anthocyanin) and the recessive sh2 allele (supersweet) lie about 0.1 Mb
apart on chr3 (0.07-0.09 cM), so a purple supersweet line needs a rare crossover inside this
interval (Anirban et al. 2023 Sci Rep 13:1050, R0002; Ravichandran et al. 2026 Plant Breed,
R0736). This step measures the interval in every assembled genome, checks gene order and
orientation, places the one published crossover, and computes how many kernels a screen needs.

Inputs
  data/raw/pan-zea.v4.pan-genes_table.tsv.gz          orthologues of the four B73 v5 genes
  data/raw/gff/<prefix>.gff3.gz                        gene coordinates and strands per genome
  results/tables/kmer_hits_near_locus.csv.gz           sh2 CDS 31-mer anchors, core genomes (step 03)
  results/tables/parents_kmer_hits_near_locus.csv.gz   same for the 8 CAAS founders (step 10)
  data/curated/genomes.csv, data/curated/hybrid_parents.csv
  data/curated/R0736_TableS1_F2.csv                    F2 kernel-colour counts, Ravichandran 2026 Table S1
  data/raw/genomes/<assembly>.fa.gz                    N runs (assembly gaps) inside the interval
                                                       (override the folder with GENOME_DIR)

Outputs
  results/tables/a1_sh2_interval_genes.csv       one row per genome x gene (sh2 = CDS anchor blocks)
  results/tables/a1_sh2_interval_summary.csv     one row per genome
  results/tables/a1_sh2_interval_gaps.csv        assembly gaps (N runs >= 10 bp) in each window
  results/tables/a1_sh2_screening_sizes.csv      kernels to screen for one purple shrunken kernel
  results/tables/a1_sh2_R0736_F2_check.csv       observed purple sh2 F2 plants in R0736 vs expectation
  results/figures/Fig9_a1_sh2_interval.png
"""
import gzip
import math
import os
import re
from multiprocessing import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon

import style

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
CUR = ROOT / "data" / "curated"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
GENOME_DIR = Path(os.environ.get("GENOME_DIR", RAW / "genomes"))
K = 31
FLANK = 50_000          # gap scan reaches 50 kb beyond a1 and sh2

GENES = {"a1": "Zm00001eb159020", "yz1": "Zm00001eb159030", "x1": "Zm00001eb159040",
         "sh2": "Zm00001eb159060"}
# Crossover in the purple-pericarp supersweet line Tim1, B73 NAM-5.0 coordinates (Anirban 2023, Fig. 3d)
TIM1_XO = (221_792_700, 221_792_792)
# a1-sh2 genetic distances reported in aleurone- and pericarp-pigmented crosses (cM)
CM_REPORTED = {"Yao et al. 2002 (PNAS)": 0.07, "Anirban et al. 2023": 0.08, "Civardi et al. 1994": 0.09}

C_A1 = style.CATEGORIES["3_pigment"][1]          # green: anthocyanin gene
C_SH2 = style.CATEGORIES["1_starch_sugar"][1]    # blue: starch gene, collinear with B73
C_INV = style.STRAND_OPP                         # orange: inverted relative to B73
C_OTHER = style.DATA_GREY_LIGHT
GROUP_LABEL = {"B73": "Temperate dent, SS", "Mo17": "Temperate dent, Lancaster", "SK": "Tropical",
               "CML247": "Tropical (CIMMYT)", "Ia453-sh2": "Sweet, supersweet", "Il14H": "Sweet, sugary",
               "P39": "Sweet, sugary", "HP301": "Popcorn", "K0326Y": "Quality-protein maize"}


def genome_table():
    core = pd.read_csv(CUR / "genomes.csv")
    core["group"] = core.corn_type + " (" + core.kernel_or_group.str.replace(r";.*", "", regex=True) + ")"
    par = pd.read_csv(CUR / "hybrid_parents.csv")
    par = par[par.in_core_set == "no"].copy()
    par["group"] = "Chinese founder, " + par.heterotic_group.str.replace(r" \(.*", "", regex=True)
    cols = ["line", "prefix", "assembly", "group"]
    return pd.concat([core[cols], par[cols]], ignore_index=True)


def orthologues(prefixes):
    """B73 symbol -> {prefix: [transcript IDs]} from the Pan-Zea v4 pan-gene table."""
    out = {s: {} for s in GENES}
    want = {v: s for s, v in GENES.items()}
    with gzip.open(RAW / "pan-zea.v4.pan-genes_table.tsv.gz", "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        col = {p: head.index(p) for p in prefixes}
        b73 = head.index("Zm00001eb")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ids = {re.sub(r"_T\d+$", "", x) for x in f[b73].split(",") if x}
            for gid in ids & set(want):
                for p, j in col.items():
                    out[want[gid]][p] = [x for x in f[j].split(",") if x and x != "NONE"]
    return out


def gff_cds(prefix, transcripts):
    """gene ID -> (chrom, CDS start, CDS end, strand) for the Pan-Zea transcript of each gene.
    CDS extents, not gene-model ends, so that UTR length and read-through models (e.g. SK x1)
    do not differ between annotations; falls back to the gene feature if no CDS is found."""
    tx = set(transcripts)
    genes = {re.sub(r"_T\d+$", "", t) for t in tx}
    gene_rows, cds = {}, {}
    with gzip.open(RAW / "gff" / f"{prefix}.gff3.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            if f[2] == "gene":
                m = re.search(r"ID=(?:gene:)?([^;\s]+)", f[8])
                if m and m.group(1) in genes:
                    gene_rows[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
            elif f[2] == "CDS":
                m = re.search(r"Parent=(?:transcript:)?([^;\s]+)", f[8])
                if m and m.group(1) in tx:
                    c = cds.setdefault(m.group(1), [f[0], int(f[3]), int(f[4]), f[6]])
                    c[1], c[2] = min(c[1], int(f[3])), max(c[2], int(f[4]))
    out = {}
    for t in tx:
        g = re.sub(r"_T\d+$", "", t)
        if t in cds:
            out[g] = tuple(cds[t])
        elif g in gene_rows:
            out[g] = gene_rows[g]
    return out


def window_gaps(job):
    """N runs >= 10 bp in chrom:lo-hi (1-based) of one assembly, streamed from the .fa.gz."""
    line, assembly, chrom, lo, hi = job
    keep, pos, buf = False, 0, []
    with gzip.open(GENOME_DIR / f"{assembly}.fa.gz", "rt") as fh:
        for row in fh:
            if row[0] == ">":
                if keep:
                    break
                keep, pos = row[1:].split()[0] == chrom, 0
                continue
            if not keep:
                continue
            seq = row.rstrip()
            a, b = pos + 1, pos + len(seq)
            if b >= lo and a <= hi:
                buf.append(seq[max(0, lo - a):hi - a + 1])
            pos = b
            if pos >= hi:
                break
    win = "".join(buf)
    gaps = [(lo + m.start(), lo + m.end() - 1) for m in re.finditer(r"[Nn]{10,}", win)]
    return line, len(win), gaps


def anchor_blocks(h):
    """Runs of sh2 CDS anchors that keep orientation and lie within 20 kb (as in step 08)."""
    if h.empty:
        return []
    h = h.sort_values("offset").assign(same=lambda d: d.qfwd == d.gfwd)
    blocks, cur = [], [h.iloc[0]]
    for _, r in h.iloc[1:].iterrows():
        p = cur[-1]
        if r.same == p.same and abs(r.pos - p.pos) < 20_000:
            cur.append(r)
        else:
            blocks.append(cur)
            cur = [r]
    blocks.append(cur)
    out = []
    for b in blocks:
        if len(b) < 10:
            continue
        d = pd.DataFrame(b)
        out.append(dict(chrom=d.chrom.iloc[0], start=int(d.pos.min()), end=int(d.pos.max()) + K - 1,
                        cds_from=int(d.offset.min()), cds_to=int(d.offset.max()) + K - 1,
                        cds_strand="+" if d.same.iloc[0] else "-"))
    return out


def screening_sizes():
    """F1 = A1 Sh2 / a1 sh2 selfed. A kernel is purple (A1 in any endosperm copy) and shrunken
    (sh2 in all three copies) only if both gametes carry sh2 and at least one carries A1:
    p = (r/2)^2 + 2 (r/2)(1 - r)/2 = (r/2)(1 - r/2). Assumes every other aleurone gene is functional
    and fixed (C1, R1, ...)."""
    rows = []
    for src, cm in CM_REPORTED.items():
        r = cm / 100
        p = (r / 2) * (1 - r / 2)
        rows.append(dict(source=src, a1_sh2_cM=cm, purple_shrunken_kernel_freq=p,
                         one_in_kernels=round(1 / p),
                         kernels_for_95pct=math.ceil(math.log(0.05) / math.log(1 - p)),
                         kernels_for_99pct=math.ceil(math.log(0.01) / math.log(1 - p))))
    return pd.DataFrame(rows)


def r0736_check(cm=0.08, kernels_per_ear=(200, 400, 600)):
    """Purple, sh2-fixed plants among the 1,843 F2 plants of R0736 (Table S1) against two expectations
    at r = 0.08 cM, with every other aleurone gene assumed functional:
      inherited  - F2 plants carrying an A1 sh2 chromosome from the F1: N * 2 * (r/2) = N * r
      ear screen - heterozygous A1 Sh2 / a1 sh2 plants (~N/2) whose selfed ear of k kernels shows at
                   least one purple shrunken kernel: N/2 * (1 - (1 - p)^k), p = (r/2)(1 - r/2)"""
    t = pd.read_csv(CUR / "R0736_TableS1_F2.csv")
    obs = int(t.filter(like="purple_").to_numpy().sum())
    n = int(t.f2_plants.sum())
    r = cm / 100
    p = (r / 2) * (1 - r / 2)
    rows = [dict(scenario="inherited recombinant chromosome", kernels_per_ear=np.nan, f2_plants=n,
                 observed_purple_sh2_plants=obs, expected=round(n * r, 2))]
    for k in kernels_per_ear:
        rows.append(dict(scenario="purple shrunken kernel on the selfed F2 ear", kernels_per_ear=k,
                         f2_plants=n, observed_purple_sh2_plants=obs,
                         expected=round(n / 2 * (1 - (1 - p) ** k), 1)))
    return pd.DataFrame(rows)


def build():
    gen = genome_table()
    orth = orthologues(gen.prefix.tolist())
    hits = pd.read_csv(TAB / "kmer_hits_near_locus.csv.gz")
    phit_path = TAB / "parents_kmer_hits_near_locus.csv.gz"
    if phit_path.exists():
        hits = pd.concat([hits, pd.read_csv(phit_path)], ignore_index=True)
    else:
        print("note: parents_kmer_hits_near_locus.csv.gz missing (run step 10); "
              "sh2 in the CAAS founders falls back to the annotated CDS")
    hits = hits[hits.v5_id == GENES["sh2"]]

    genes, summ = [], []
    for g in gen.itertuples():
        tx = [t for s in GENES for t in orth[s].get(g.prefix, [])]
        coords = gff_cds(g.prefix, tx)
        loc = {}
        for s in ("a1", "yz1", "x1", "sh2"):
            ids = [re.sub(r"_T\d+$", "", t) for t in orth[s].get(g.prefix, [])]
            cand = [coords[i] + (i,) for i in ids if i in coords]
            if s != "a1" and "a1" in loc:          # the copy next to a1 if there are several
                cand = sorted([c for c in cand if c[0] == loc["a1"][0]],
                              key=lambda c: abs(c[1] - loc["a1"][1]))
            if cand:
                loc[s] = cand[0]
        if "a1" not in loc or "sh2" not in loc:
            print(f"{g.line}: a1 or sh2 orthologue not found")
            continue
        chrom = loc["a1"][0]
        a1s, a1e = loc["a1"][1], loc["a1"][2]
        sh2_mid = (loc["sh2"][1] + loc["sh2"][2]) / 2
        sign = 1 if sh2_mid > a1e else -1
        x0 = a1e if sign == 1 else a1s                     # a1 CDS end facing sh2 = 0 kb

        def rel(p):
            return (p - x0) * sign

        def ostrand(st):                                   # strand in the oriented frame
            return st if sign == 1 else {"+": "-", "-": "+"}[st]

        for s in ("a1", "yz1", "x1"):
            if s in loc:
                c, st, en, strand, gid = loc[s]
                genes.append(dict(line=g.line, symbol=s, gene_id=gid, chrom=c, start=st, end=en,
                                  strand=strand, x_from_kb=min(rel(st), rel(en)) / 1e3,
                                  x_to_kb=max(rel(st), rel(en)) / 1e3, oriented_strand=ostrand(strand),
                                  source="CDS of the Pan-Zea v4 orthologue transcript"))
        blocks = anchor_blocks(hits[hits.line == g.line])
        if blocks:
            src = "sh2 CDS 31-mer anchors"
        else:                                              # fall back to the annotated CDS
            c, st, en, strand, gid = loc["sh2"]
            blocks = [dict(chrom=c, start=st, end=en, cds_from=np.nan, cds_to=np.nan,
                           cds_strand="-" if strand == "-" else "+")]
            src = "annotated CDS (no anchors)"
        b73_like = []
        for b in blocks:
            # B73 sh2 is on the minus strand: an anchor block on the opposite strand to the CDS
            # in a genome where a1 precedes sh2 runs the same way as in B73
            genome_strand = "-" if b["cds_strand"] == "-" else "+"
            collinear = ostrand(genome_strand) == "-"
            b73_like.append(collinear)
            genes.append(dict(line=g.line, symbol="sh2", gene_id=loc["sh2"][4], chrom=b["chrom"],
                              start=b["start"], end=b["end"], strand=genome_strand,
                              x_from_kb=min(rel(b["start"]), rel(b["end"])) / 1e3,
                              x_to_kb=max(rel(b["start"]), rel(b["end"])) / 1e3,
                              oriented_strand=ostrand(genome_strand), cds_from=b["cds_from"],
                              cds_to=b["cds_to"], collinear_with_B73=collinear, source=src))
        gx = pd.DataFrame([r for r in genes if r["line"] == g.line])
        sh = gx[gx.symbol == "sh2"]
        others = gx[gx.symbol.isin(["yz1", "x1"])]
        summ.append(dict(
            line=g.line, group=g.group, chrom=chrom, a1_start=a1s, a1_end=a1e,
            genes_in_order="-".join(gx.sort_values("x_from_kb").symbol.drop_duplicates()),
            yz1_x1_same_strand_as_a1=bool((others.strand == loc["a1"][3]).all()) if len(others) else None,
            a1_to_yz1_kb=round(gx.loc[gx.symbol == "yz1", "x_from_kb"].min(), 1) if "yz1" in loc else None,
            a1_to_x1_kb=round(gx.loc[gx.symbol == "x1", "x_from_kb"].min(), 1) if "x1" in loc else None,
            a1_to_sh2_cds_kb=round(sh.x_from_kb.min(), 1),
            a1_to_sh2_far_end_kb=round(sh.x_to_kb.max(), 1),
            sh2_blocks=len(sh),
            sh2_structure=("split; part inverted relative to B73" if (len(sh) > 1 and not all(b73_like))
                           else "collinear with B73" if all(b73_like) else "inverted relative to B73"),
            sh2_source=src))
    genes, summ = pd.DataFrame(genes), pd.DataFrame(summ)

    # assembly gaps (N runs): is the a1-sh2 distance read from contiguous sequence?
    asm = gen.set_index("line").assembly
    jobs = [(r.line, asm[r.line], r.chrom,
             max(1, int(genes.loc[genes.line == r.line, "start"].min()) - FLANK),
             int(genes.loc[genes.line == r.line, "end"].max()) + FLANK) for r in summ.itertuples()]
    with Pool(int(os.environ.get("GAP_WORKERS", 4))) as pool:
        res = {line: (n, gaps) for line, n, gaps in pool.map(window_gaps, jobs)}
    extra, gap_rows = [], []
    for r, job in zip(summ.itertuples(), jobs):
        n, gaps = res[r.line]
        assert n == job[4] - job[3] + 1, f"{r.line}: window not read in full"
        gx = genes[genes.line == r.line]
        sh = gx[gx.symbol == "sh2"]
        lo, hi = min(r.a1_start, sh.start.min()), max(r.a1_end, sh.end.max())
        inside = [g for g in gaps if g[1] >= lo and g[0] <= hi]
        sign = 1 if sh.start.mean() > r.a1_end else -1
        x0 = r.a1_end if sign == 1 else r.a1_start
        for a, b in gaps:
            xs = sorted([(a - x0) * sign / 1e3, (b - x0) * sign / 1e3])
            gap_rows.append(dict(line=r.line, chrom=r.chrom, start=a, end=b, length_bp=b - a + 1,
                                 x_from_kb=xs[0], x_to_kb=xs[1], inside_a1_sh2=(b >= lo and a <= hi)))
        near = sh.iloc[(sh.x_from_kb).argmin()]                # sh2 block nearest to a1
        left = max([g[1] for g in gaps if g[1] < near.start], default=None)
        right = min([g[0] for g in gaps if g[0] > near.end], default=None)
        seg_lo = left + 1 if left else job[3]
        seg_hi = right - 1 if right else job[4]
        extra.append(dict(interval_gaps=len(inside), interval_gap_bp=sum(b - a + 1 for a, b in inside),
                          sh2_gapfree_segment_kb=round((seg_hi - seg_lo + 1) / 1e3, 1),
                          sh2_segment_bounded_by_gaps=bool(left and right)))
    summ = pd.concat([summ, pd.DataFrame(extra)], axis=1)
    gap_df = pd.DataFrame(gap_rows, columns=["line", "chrom", "start", "end", "length_bp", "x_from_kb",
                                             "x_to_kb", "inside_a1_sh2"])
    summ["assembly_support"] = np.where(summ.interval_gaps == 0, "contiguous sequence",
                                        "interval contains assembly gaps; distance depends on scaffolding")

    # the published crossover, in B73 v5 coordinates
    b = summ.set_index("line").loc["B73"]
    xo_kb = ((TIM1_XO[0] - b.a1_end) / 1e3, (TIM1_XO[1] - b.a1_end) / 1e3)
    yz1 = genes[(genes.line == "B73") & (genes.symbol == "yz1")].iloc[0]
    print(f"Tim1 crossover (Anirban 2023): B73 chr3:{TIM1_XO[0]:,}-{TIM1_XO[1]:,} = "
          f"{xo_kb[0]:.1f}-{xo_kb[1]:.1f} kb distal of a1, {(yz1.start - TIM1_XO[1]) / 1e3:.1f} kb "
          f"before yz1 ({100 * xo_kb[0] / b.a1_to_sh2_cds_kb:.0f}% of the a1-sh2 gap)")
    return genes, summ, gap_df, xo_kb


def gene_arrow(ax, x0, x1, y, strand, color, h=0.42, hatch=None, edge=None):
    head = min(0.45 * (x1 - x0), 2.2)
    if strand == "+":
        pts = [(x0, y - h / 2), (x1 - head, y - h / 2), (x1, y), (x1 - head, y + h / 2), (x0, y + h / 2)]
    else:
        pts = [(x1, y - h / 2), (x0 + head, y - h / 2), (x0, y), (x0 + head, y + h / 2), (x1, y + h / 2)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor=edge or color, lw=0.4,
                         hatch=hatch, zorder=3))


def figure(genes, summ, gap_df, xo_kb):
    style.setup()
    order = ["B73", "Mo17", "SK", "CML247", "Ia453-sh2", "Il14H", "P39", "HP301", "K0326Y",
             "Huangzaosi", "Chang7-2", "Jing92", "Dan340", "Ye478", "Zheng58", "Jing724", "Xu178"]
    order = [l for l in order if l in set(summ.line)] + sorted(set(summ.line) - set(order))
    n = len(order)
    H_MM = 112
    fig = plt.figure(figsize=(style.WIDTH_IN, H_MM / style.MM))
    left, right = 1.38, 0.95                                        # inches for label columns
    ax = style.axes_in(fig, left, 0.62, style.WIDTH_IN - left - right, H_MM / style.MM - 0.62 - 0.30)
    XCAP = 180                     # kb; a gene placed further away is drawn after an axis break
    xmax = XCAP + 12
    ax.set_xlim(-6, xmax)
    ax.set_ylim(n - 0.4, -1.25)
    style.despine(ax, keep=("bottom",))
    ax.set_yticks([])
    ax.set_xlabel("Distance from " + style.italic_gene("a1") + " CDS (kb)")
    ax.set_xticks(range(0, XCAP + 1, 20))
    ax.set_xticks(range(0, XCAP + 1, 10), minor=True)
    ax.spines["bottom"].set_bounds(-6, XCAP)

    S = summ.set_index("line")
    trans_lab = ax.get_yaxis_transform()
    for i, line in enumerate(order):
        gx = genes[genes.line == line].copy()
        far = gx.x_from_kb > XCAP
        if far.any():                                      # axis break: draw the far gene at the edge
            shift = gx.loc[far, "x_from_kb"].min() - (XCAP + 4)
            far_kb = gx.loc[far, "x_from_kb"].min()
            gx.loc[far, ["x_from_kb", "x_to_kb"]] -= shift
            ax.plot([gx.x_from_kb.min() - 1, XCAP - 1], [i, i], color=style.CHR_EDGE, lw=0.5, zorder=1)
            ax.plot([XCAP + 3, gx.x_to_kb.max() + 1], [i, i], color=style.CHR_EDGE, lw=0.5, zorder=1)
            for dx in (-0.6, 1.4):
                ax.plot([XCAP + dx - 0.7, XCAP + dx + 0.7], [i + 0.3, i - 0.3], color=style.INK2, lw=0.5)
            ax.text(XCAP + 1, i - 0.36, f"{far_kb / 1e3:.2f} Mb", fontsize=style.FS_SMALL, ha="center",
                    va="bottom", color=style.INK2)
        else:
            ax.plot([gx.x_from_kb.min() - 1, gx.x_to_kb.max() + 1], [i, i], color=style.CHR_EDGE,
                    lw=0.5, zorder=1, solid_capstyle="butt")
        # assembly gaps in the drawn range: short black bars across the track
        gp = gap_df[(gap_df.line == line) & (gap_df.x_to_kb > -6) & (gap_df.x_from_kb < XCAP)]
        for gr in gp.itertuples():
            ax.plot([max(gr.x_from_kb, -6), min(gr.x_to_kb, XCAP)], [i, i], color=style.INK, lw=1.6,
                    solid_capstyle="butt", zorder=2)
        for r in gx.itertuples():
            if r.symbol == "a1":
                c, hatch = C_A1, None
            elif r.symbol == "sh2":
                c, hatch = (C_SH2, None) if r.collinear_with_B73 else (C_INV, None)
            else:
                c, hatch = C_OTHER, None
            gene_arrow(ax, r.x_from_kb, r.x_to_kb, i, r.oriented_strand, c, hatch=hatch)
        sh = gx[gx.symbol == "sh2"].sort_values("x_from_kb")
        if len(sh) > 1:                                    # split gene: join the parts
            xa, xb = sh.x_to_kb.iloc[0], sh.x_from_kb.iloc[1]
            xs = np.linspace(xa, xb, 40)
            ax.plot(xs, i - 0.22 - 0.28 * np.sin(np.pi * (xs - xa) / (xb - xa)), color=style.INK2,
                    lw=0.5, zorder=2)
            ax.text((xa + xb) / 2, i - 0.55, f"{xb - xa:.1f} kb", fontsize=style.FS_SMALL,
                    ha="center", va="bottom", color=style.INK2)
        # left labels: line and group; right label: a1 -> sh2 gap
        bold = "bold" if line in ("B73", "Ia453-sh2") else "normal"
        ax.text(-0.012, i, line, transform=trans_lab, ha="right", va="center", fontsize=style.FS,
                fontweight=bold)
        grp = GROUP_LABEL.get(line, S.loc[line, "group"])
        ax.text(-0.012, i + 0.36, grp, transform=trans_lab, ha="right", va="center",
                fontsize=style.FS_SMALL - 0.5, color=style.INK2)
        gap_kb = S.loc[line, "a1_to_sh2_cds_kb"]
        mark = ("†" if S.loc[line, "interval_gaps"] else "") + \
               ("*" if S.loc[line, "sh2_source"] != "sh2 CDS 31-mer anchors" else "")
        ax.text(1.012, i, (f"{gap_kb:,.1f}" if gap_kb < 1000 else f"{gap_kb / 1e3:.2f} Mb") + mark,
                transform=trans_lab, ha="left", va="center", fontsize=style.FS)
    ax.text(1.012, -0.95, style.italic_gene("a1–sh2") + "\ngap (kb)", transform=trans_lab, ha="left",
            va="center", fontsize=style.FS, linespacing=1.0)

    # gene names over the reference track
    b = genes[genes.line == order[0]]
    for r in b.drop_duplicates("symbol").itertuples():
        x = (r.x_from_kb + b[b.symbol == r.symbol].x_to_kb.max()) / 2
        ax.text(x, -0.62, style.italic_gene(r.symbol), ha="center", va="bottom", fontsize=style.FS)
    # Tim1 crossover on B73
    xm = sum(xo_kb) / 2
    ax.plot([xm], [-0.36], marker="v", ms=3.2, color=style.INK, zorder=4, clip_on=False)
    ax.annotate("Tim1 crossover", xy=(xm, -0.42), xytext=(xm + 9, -0.98), fontsize=style.FS_SMALL,
                ha="left", va="center", color=style.INK,
                arrowprops=dict(arrowstyle="-", lw=0.4, color=style.INK2, shrinkA=0.5, shrinkB=1.5))
    # separators between germplasm blocks
    for yb in (3.5, 8.5):
        if yb < n:
            ax.axhline(yb, color=style.CHR_FILL, lw=0.6, zorder=0)

    handles = [Patch(facecolor=C_A1, edgecolor=C_A1, label=style.italic_gene("a1") + " (anthocyaninless1)"),
               Patch(facecolor=C_OTHER, edgecolor=C_OTHER, label=style.italic_gene("yz1") + ", "
                     + style.italic_gene("x1")),
               Patch(facecolor=C_SH2, edgecolor=C_SH2, label=style.italic_gene("sh2") + " CDS, collinear with B73"),
               Patch(facecolor=C_INV, edgecolor=C_INV, label=style.italic_gene("sh2") + " CDS, inverted"),
               Line2D([], [], color=style.INK, lw=1.6, label="Assembly gap (N)"),
               Line2D([], [], marker="v", ls="", ms=3.2, color=style.INK, label="Crossover in Tim1 (B73 frame)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 0.0),
               ncol=6, fontsize=style.FS, handlelength=1.2, columnspacing=0.9, frameon=False)
    return fig


def main():
    genes, summ, gap_df, xo_kb = build()
    TAB.mkdir(parents=True, exist_ok=True)
    genes.to_csv(TAB / "a1_sh2_interval_genes.csv", index=False)
    summ.to_csv(TAB / "a1_sh2_interval_summary.csv", index=False)
    gap_df.to_csv(TAB / "a1_sh2_interval_gaps.csv", index=False)
    scr = screening_sizes()
    scr.to_csv(TAB / "a1_sh2_screening_sizes.csv", index=False)
    chk = r0736_check()
    chk.to_csv(TAB / "a1_sh2_R0736_F2_check.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ.drop(columns=["a1_start", "a1_end"]).to_string(index=False))
        print(scr.to_string(index=False))
        print(chk.to_string(index=False))
    fig = figure(genes, summ, gap_df, xo_kb)
    FIG.mkdir(parents=True, exist_ok=True)
    style.save_png(fig, FIG / "Fig9_a1_sh2_interval.png")


if __name__ == "__main__":
    main()
