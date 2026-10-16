# -*- coding: utf-8 -*-
"""Sweet corn numbers quoted in the review, computed from the result tables in ../results/tables/.

1  Chromosome 4 diversity contrast (Hu et al. 2021 GBS calls, 1-Mb windows on B73 RefGen_v4,
   table popgen_windows_Hu2021_1Mb.csv written by 07_population_differentiation.py).
   For each window the reduction in expected heterozygosity of sweet corn relative to temperate
   field corn is 1 - He_sweet / He_temperate_field. Reported: the median and interquartile range
   (25th-75th percentile) over
     a. the windows from 40 to 160 Mb of chromosome 4 (start >= 40 Mb and end <= 160 Mb);
     b. the windows whose centre lies within 20 Mb of su1 (gene start on B73 v4 from
        popgen_key_genes.csv), with 10 and 40 Mb as a sensitivity check;
     c. all windows of the genome.
   Because the 40-160 Mb window was chosen after looking at the scan, the script also ranks it
   among all blocks of 120 consecutive 1-Mb windows (blocks that contain a window without data are
   skipped).
2  The a1-sh2 interval (a1_sh2_interval_summary.csv, written by 13_a1_sh2_interval.py): distance
   from the a1 CDS to the sh2 CDS over the assemblies in which sh2 is collinear with B73 and the
   interval holds less than 1 kb of assembly gap.
3  Chromosomes of the endosperm genes (key_gene_B73v5.csv, B73 RefGen_v5).

Writes the same text to ../results/genome_numbers.txt.
Usage: python -I -B genome/scripts/genome_numbers.py
"""
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.normpath(os.path.join(HERE, "..", "results", "tables"))
OUTFILE = os.path.normpath(os.path.join(HERE, "..", "results", "genome_numbers.txt"))
BLOCK = 120          # windows in the 40-160 Mb block
MAX_GAP_BP = 1000    # assemblies with more assembly gap inside the a1-sh2 interval are left out


def summary(x):
    return x.median(), x.quantile(0.25), x.quantile(0.75), len(x)


def pct(v):
    return "%.1f%%" % (100 * v)


def main():
    out = []
    say = out.append

    w = pd.read_csv(os.path.join(TAB, "popgen_windows_Hu2021_1Mb.csv"))
    w["reduction"] = 1 - w.he_sweet / w.he_temperate_field
    w["mid"] = (w.start + w.end) / 2
    kg = pd.read_csv(os.path.join(TAB, "popgen_key_genes.csv"))
    su1 = kg.loc[kg.symbol == "su1"].iloc[0]
    su1_pos, su1_chr = int(su1.v4_start), int(su1.chr)
    c4 = w[w.chr == su1_chr]

    say("Genome numbers for the review (genome_numbers.py)")
    say("")
    say("1. Chromosome 4: reduction in expected heterozygosity of sweet corn relative to temperate field "
        "corn, 1 - He_sweet / He_field, per 1-Mb window (Hu et al. 2021 GBS, B73 v4)")
    a = c4[(c4.start >= 40e6) & (c4.end <= 160e6)].reduction
    m, q1, q3, n = summary(a)
    say("   a. chromosome 4, 40-160 Mb: median %s (IQR %s to %s), %d windows" % (pct(m), pct(q1), pct(q3), n))
    for r in (20, 10, 40):
        b = c4[(c4.mid - su1_pos).abs() <= r * 1e6].reduction
        m, q1, q3, n = summary(b)
        tag = "b." if r == 20 else "  "
        say("   %s windows centred within %d Mb of su1 (chr%d:%s, B73 v4): median %s (IQR %s to %s), %d windows"
            % (tag, r, su1_chr, format(su1_pos, ","), pct(m), pct(q1), pct(q3), n))
    m, q1, q3, n = summary(w.reduction)
    say("   c. all windows of the genome: median %s (IQR %s to %s), %s windows"
        % (pct(m), pct(q1), pct(q3), format(n, ",")))

    rows = []
    for c, g in w.groupby("chr"):
        red = g.set_index("win").reduction
        have = set(red.index)
        for k in range(int(red.index.min()), int(red.index.max()) - BLOCK + 2):
            idx = list(range(k, k + BLOCK))
            if all(i in have for i in idx):
                rows.append((c, k, red.loc[idx].median()))
    blocks = pd.DataFrame(rows, columns=["chr", "start_mb", "median"])
    best = blocks["median"].max()
    top = blocks[blocks["median"] >= best - 1e-12]
    own = blocks[(blocks.chr == su1_chr) & (blocks.start_mb == 40)]["median"].iloc[0]
    other = blocks[blocks.chr != su1_chr].sort_values("median").iloc[-1]
    missing = [(int(c), sorted(set(range(int(g.win.min()), int(g.win.max()) + 1)) - set(g.win)))
               for c, g in w.groupby("chr")]
    missing = [(c, m_) for c, m_ in missing if m_]
    say("   Rank of the 40-160 Mb block among the %d blocks of %d consecutive 1-Mb windows (blocks with a "
        "window lacking data skipped; windows without data: %s): its median %s is the largest in the genome "
        "(shared by blocks starting at %s Mb on chromosome %d); the largest on any other chromosome is %s "
        "(chromosome %d, %d-%d Mb)."
        % (len(blocks), BLOCK, "; ".join("chr%d window %s" % (c, ",".join(str(x) for x in m_)) for c, m_ in missing),
           pct(own), ", ".join(str(int(x)) for x in sorted(top[top.chr == su1_chr].start_mb)), su1_chr,
           pct(other["median"]), other.chr, other.start_mb, other.start_mb + BLOCK))
    assert abs(own - best) < 1e-12

    # 2 ---------------------------------------------------------------------------------------------
    s = pd.read_csv(os.path.join(TAB, "a1_sh2_interval_summary.csv"))
    keep = s[(s.sh2_structure == "collinear with B73") & (s.interval_gap_bp < MAX_GAP_BP)]
    drop = s[~s.index.isin(keep.index)]
    lo, hi = keep.loc[keep.a1_to_sh2_cds_kb.idxmin()], keep.loc[keep.a1_to_sh2_cds_kb.idxmax()]
    say("")
    say("2. a1-sh2 interval (a1 CDS end to sh2 CDS), assemblies with sh2 collinear with B73 and < 1 kb of "
        "assembly gap in the interval")
    say("   %d assemblies: %.1f kb (%s) to %.1f kb (%s); median %.1f kb"
        % (len(keep), lo.a1_to_sh2_cds_kb, lo.line, hi.a1_to_sh2_cds_kb, hi.line, keep.a1_to_sh2_cds_kb.median()))
    say("   included: %s" % ", ".join(keep.line))
    for r in drop.itertuples():
        say("   left out: %s (%s; %d assembly gaps, %s bp; distance %.1f kb)"
            % (r.line, r.sh2_structure, r.interval_gaps, format(int(r.interval_gap_bp), ","), r.a1_to_sh2_cds_kb))
    gapped = keep[keep.interval_gaps > 0]
    for r in gapped.itertuples():
        say("   kept with a small gap: %s (%d gap, %d bp); without it the range is unchanged: %s"
            % (r.line, r.interval_gaps, r.interval_gap_bp,
               "yes" if lo.line != r.line and hi.line != r.line else "no"))

    # 3 ---------------------------------------------------------------------------------------------
    g = pd.read_csv(os.path.join(TAB, "key_gene_B73v5.csv"))
    say("")
    say("3. Chromosomes of endosperm and linked genes (B73 RefGen_v5)")
    for sym in ("a1", "sh2", "bt2", "su1", "se1", "wx1"):
        r = g[g.symbol == sym].iloc[0]
        say("   %-4s %s:%s-%s (%.1f Mb)" % (sym, r.chrom, format(int(r.start), ","), format(int(r.end), ","),
                                          r.start / 1e6))

    text = "\n".join(out) + "\n"
    print(text, end="")
    with open(OUTFILE, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print("written to %s" % OUTFILE)


if __name__ == "__main__":
    main()
