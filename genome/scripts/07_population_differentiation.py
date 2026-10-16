"""Genome-wide footprints of food-corn breeding: Fst and diversity scans.

Data (both CC BY; copied into data/raw/popgen by script 01):
  Hu et al. 2021 - GBS VCF, 5,643 samples x ~860k SNPs on B73 v4, with the accession groups
                   of Supplementary Data 4 (798 sweet corn, 265 stiff stalk, 255 non-stiff stalk,
                   566 tropical, 126 popcorn ...). Sample names are matched by replacing ':' with
                   '_' and, if needed, dropping the leading project tag (e.g. 'AD_').
  Qu et al. 2026 - 5,759-SNP GBTS array (Table S4): 51 sweet, 72 waxy, 75 field lines, B73 v4.

Statistics
  Hudson Fst per SNP, combined over windows as a ratio of averages (Bhatia et al. 2013):
    num = (p1-p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1);  den = p1(1-p2) + p2(1-p1)
  Expected heterozygosity He = 2p(1-p) n/(n-1), averaged over SNPs in a window.
  SNP filters: >= MIN_ALLELES called alleles in every group of a contrast; pooled MAF >= MIN_MAF.

Outputs
  results/intermediate/Hu2021_group_allele_counts.csv.gz   per-SNP counts per group
  results/tables/popgen_windows_Hu2021_1Mb.csv              window Fst / He
  results/tables/popgen_windows_Qu2026_5Mb.csv / popgen_snps_Qu2026.csv (per-SNP Fst)
  results/tables/popgen_key_genes.csv                       statistics at each key gene + percentile
  results/tables/popgen_top_windows_sweet_vs_field.csv      top 1% windows with the genes inside
  results/figures/Fig5_genome_scans_food_corn_selection.png (600 dpi, 183 mm wide)
"""
import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import style

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
POP = RAW / "popgen"
TAB = ROOT / "results" / "tables"
INT = ROOT / "results" / "intermediate"
FIG = ROOT / "results" / "figures"
INT.mkdir(parents=True, exist_ok=True)

HU_GROUPS = {  # analysis group -> Supplementary Data 4 GroupID values
    "sweet": ["Sweet Corn"],
    "temperate_field": ["stiff stalk", "non-stiff stalk"],
    "tropical": ["tropical"],
    "popcorn": ["popcorn"],
}
HU_CONTRASTS = [("sweet", "temperate_field"), ("popcorn", "temperate_field"), ("sweet", "tropical")]
QU_CONTRASTS = [("Waxy", "Field"), ("Sweet", "Field"), ("Sweet", "Waxy")]
MIN_ALLELES = 40
MIN_MAF = 0.02
WIN_HU = 1_000_000
WIN_QU = 5_000_000
QU_FLANK = 2_000_000   # per-SNP search radius around a key gene (array density ~1 SNP / 370 kb)
BATCH = 2000


# ---- helpers ----------------------------------------------------------------------
def hudson(p1, n1, p2, n2):
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return num, den


def window_stats(df, groups, contrasts, win):
    """df: chr, pos, alt_<g>, n_<g>. Returns per-window Fst (ratio of averages) and mean He."""
    df = df.copy()
    df["win"] = (df.pos // win).astype(int)
    out = []
    for g in groups:
        p = df[f"alt_{g}"] / df[f"n_{g}"]
        df[f"he_{g}"] = 2 * p * (1 - p) * df[f"n_{g}"] / (df[f"n_{g}"] - 1)
    for a, b in contrasts:
        ok = (df[f"n_{a}"] >= MIN_ALLELES) & (df[f"n_{b}"] >= MIN_ALLELES)
        pa, pb = df[f"alt_{a}"] / df[f"n_{a}"], df[f"alt_{b}"] / df[f"n_{b}"]
        pool = (df[f"alt_{a}"] + df[f"alt_{b}"]) / (df[f"n_{a}"] + df[f"n_{b}"])
        ok &= np.minimum(pool, 1 - pool) >= MIN_MAF
        num, den = hudson(pa, df[f"n_{a}"], pb, df[f"n_{b}"])
        df[f"num_{a}_{b}"] = np.where(ok, num, np.nan)
        df[f"den_{a}_{b}"] = np.where(ok, den, np.nan)
        df[f"ok_{a}_{b}"] = ok
    agg = {"pos": "count"}
    for a, b in contrasts:
        agg[f"num_{a}_{b}"] = "sum"
        agg[f"den_{a}_{b}"] = "sum"
        agg[f"ok_{a}_{b}"] = "sum"
    for g in groups:
        agg[f"he_{g}"] = "mean"
    w = df.groupby(["chr", "win"]).agg(agg).reset_index().rename(columns={"pos": "n_snps"})
    for a, b in contrasts:
        w[f"fst_{a}_vs_{b}"] = w[f"num_{a}_{b}"] / w[f"den_{a}_{b}"]
        w.loc[w[f"ok_{a}_{b}"] < 5, f"fst_{a}_vs_{b}"] = np.nan   # too few informative SNPs
        w = w.drop(columns=[f"num_{a}_{b}", f"den_{a}_{b}"]).rename(columns={f"ok_{a}_{b}": f"nsnp_{a}_vs_{b}"})
    w["start"] = w.win * win
    w["end"] = w.start + win
    return w


# ---- Hu et al. 2021 VCF --------------------------------------------------------------
def hu_counts():
    cache = INT / "Hu2021_group_allele_counts.csv.gz"
    if cache.exists():
        return pd.read_csv(cache)
    grp = pd.read_excel(POP / "Hu2021_SupplementaryData4_groups.xlsx", sheet_name="TreeKey", header=1)
    grp = dict(zip(grp["Accession ID"].astype(str), grp["GroupID"]))
    vcf = POP / "Hu2021_noSeed_denovo_publicSweet.vcf.gz"
    with gzip.open(vcf, "rb") as fh:
        for line in fh:
            if line.startswith(b"#CHROM"):
                samples = line.rstrip(b"\n").decode().split("\t")[9:]
                break

        def norm(s):
            a = s.replace(":", "_")
            if a in grp:
                return a
            b = a.split("_", 1)[1] if "_" in a else a
            return b if b in grp else None
        label = [grp.get(norm(s)) if norm(s) else None for s in samples]
        idx = {g: np.array([i for i, l in enumerate(label) if l in v]) for g, v in HU_GROUPS.items()}
        for g, ix in idx.items():
            print(f"  Hu2021 group {g}: {len(ix)} samples")
        n = len(samples)
        width = 4 * n - 1
        rows = []
        chrs, poss, gts = [], [], []

        def flush():
            if not gts:
                return
            buf = np.frombuffer(b"\t".join(gts) + b"\t", dtype=np.uint8).reshape(len(gts), n, 4)
            a1, a2 = buf[:, :, 0], buf[:, :, 2]
            alt = (a1 > 48).astype(np.int16) + (a2 > 48).astype(np.int16)        # '1'..'9' = alt
            called = (a1 != 46).astype(np.int16) + (a2 != 46).astype(np.int16)   # '.' = missing
            d = {"chr": chrs[:], "pos": poss[:]}
            for g, ix in idx.items():
                d[f"alt_{g}"] = alt[:, ix].sum(1)
                d[f"n_{g}"] = called[:, ix].sum(1)
            rows.append(pd.DataFrame(d))
            chrs.clear(), poss.clear(), gts.clear()

        skipped = 0
        for line in fh:
            f = line.rstrip(b"\n").split(b"\t", 9)
            if len(f[9]) != width:      # non-GT or multi-character calls: skip
                skipped += 1
                continue
            chrs.append(f[0].decode())
            poss.append(int(f[1]))
            gts.append(f[9])
            if len(gts) == BATCH:
                flush()
        flush()
    counts = pd.concat(rows, ignore_index=True)
    counts = counts[counts.chr.isin([str(i) for i in range(1, 11)])]
    counts["chr"] = counts.chr.astype(int)
    print(f"  Hu2021: {len(counts):,} SNPs on chr1-10 ({skipped} malformed lines skipped)")
    counts.to_csv(cache, index=False, compression={"method": "gzip", "mtime": 0})
    return counts


# ---- Qu et al. 2026 array -------------------------------------------------------------
def qu_counts_shared_ref():
    """Per-group alt-allele counts, every group polarised against the Table S2 first allele."""
    t = pd.read_excel(POP / "Qu2026_supplementary_tables.xlsx", sheet_name="Table S4", header=2)
    t = t.rename(columns={"ID": "snp", "Chromosome": "chr", "Position": "pos"})
    t = t[pd.to_numeric(t.chr, errors="coerce").between(1, 10)].reset_index(drop=True)
    t["chr"] = t.chr.astype(int)
    t["pos"] = t.pos.astype(int)
    info = pd.read_excel(POP / "Qu2026_supplementary_tables.xlsx", sheet_name="Table S2", header=2)
    ref_allele = dict(zip(info["SNP name"], info["Alleles"].astype(str).str[0]))
    ref = t.snp.map(ref_allele).fillna("N").to_numpy()
    out = t[["snp", "chr", "pos"]].copy()
    for g in ("Sweet", "Waxy", "Field"):
        cols = [c for c in t.columns if str(c).startswith(g + "_")]
        G = t[cols].astype(str)
        c0 = G.apply(lambda col: col.str[0]).to_numpy()
        c1 = G.apply(lambda col: col.str[1]).to_numpy()
        valid = np.isin(c0, list("ACGT")) & np.isin(c1, list("ACGT"))
        r = ref[:, None]
        out[f"alt_{g}"] = (valid * ((c0 != r).astype(int) + (c1 != r).astype(int))).sum(1)
        out[f"n_{g}"] = (2 * valid).sum(1)
        print(f"  Qu2026 group {g}: {len(cols)} lines")
    return out


# ---- key genes on B73 v4 --------------------------------------------------------------
def key_genes_v4():
    genes = pd.read_csv(ROOT / "data" / "curated" / "key_genes.csv")
    x = pd.read_csv(RAW / "B73v5_to_B73v4.tsv", sep="\t", header=None, names=["v5", "v4"])
    v4map = dict(zip(x.v5, x.v4))
    fd = pd.read_csv(RAW / "Zm00001d.2.full_gene_data.txt.gz", sep="\t", dtype=str)
    fd.columns = [c.strip() for c in fd.columns]
    fd = fd.apply(lambda s: s.str.strip())
    fd = fd.set_index("gene_model")
    rows = []
    for g in genes.itertuples():
        cands = [v for v in str(v4map.get(g.v5_id, "")).split(",") if v in fd.index]
        pick = None
        for v in cands:
            c = style.chrom_number("chr" + fd.loc[v, "gene_model_chr"].replace("chr", "").replace("Chr", ""))
            if c is not None:
                pick = (v, c, int(fd.loc[v, "gene_model_start"]), int(fd.loc[v, "gene_model_end"]))
                break
        rows.append(dict(symbol=g.symbol, label=g.label, category=g.category, v5_id=g.v5_id,
                         v4_id=pick[0] if pick else None, chr=pick[1] if pick else None,
                         v4_start=pick[2] if pick else None, v4_end=pick[3] if pick else None))
    return pd.DataFrame(rows)


def v4_gene_table():
    fd = pd.read_csv(RAW / "Zm00001d.2.full_gene_data.txt.gz", sep="\t", dtype=str)
    fd.columns = [c.strip() for c in fd.columns]
    fd = fd.apply(lambda s: s.str.strip())
    fd["chr"] = pd.to_numeric(fd.gene_model_chr.str.replace("(?i)chr", "", regex=True), errors="coerce")
    fd = fd[fd.chr.between(1, 10)].copy()
    fd["start"] = fd.gene_model_start.astype(int)
    fd["end"] = fd.gene_model_end.astype(int)
    return fd[["gene_model", "chr", "start", "end", "locus_symbol", "locus_name"]]


# ---- figure ----------------------------------------------------------------------
def manhattan(ax, w, col, kg, ylabel, highlight, offsets, ylim=None, ref_line=None):
    """Scan points in alternating greys by chromosome; key genes as category-coloured lines with
    a marker (category shape) and an italic label above the axes."""
    w = w.dropna(subset=[col]).copy()
    mid = (w.start + w.end) / 2 if "start" in w else w.pos
    w["x"] = w.chr.map(offsets) + mid
    for c in range(1, 11):
        s = w[w.chr == c]
        ax.plot(s.x / 1e6, s[col], ls="none", marker="o", ms=1.1, mew=0,
                color=style.DATA_GREY if c % 2 else style.DATA_GREY_LIGHT, rasterized=True, zorder=2)
    if ref_line is None:   # F_ST panels: mark the genome-wide 99th percentile
        thr = np.nanpercentile(w[col], 99)
        ax.axhline(thr, color=style.INK2, lw=0.5, ls=(0, (3, 2)), zorder=3)
        ax.annotate("99th pct.", (1.0, thr), xycoords=("axes fraction", "data"), xytext=(2, 0),
                    textcoords="offset points", fontsize=style.FS_SMALL, va="center", color=style.INK2)
    else:
        ax.axhline(ref_line, color=style.INK2, lw=0.5, zorder=3)
    xmax = (offsets[10] + 160e6) / 1e6
    ax.set_xlim(0, xmax)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel)
    sel = kg[kg.symbol.isin(highlight)].dropna(subset=["chr"]).copy()
    sel["x"] = [(offsets[int(c)] + v) / 1e6 for c, v in zip(sel.chr, sel.v4_start)]
    sel = sel.sort_values("x")
    w_in, h_in = style.axes_size_in(ax)
    gap = 1.3 * style.FS_TICK / 72 * xmax / w_in          # Mb between neighbouring rotated labels
    x_lab = style.spread(sel.x.tolist(), gap, lo=0, hi=xmax)
    tr = ax.get_xaxis_transform()                           # x in data, y in axes fraction
    f = 1 / h_in                                            # axes fraction per inch
    for r, xl in zip(sel.itertuples(), x_lab):
        name, colr, mk = style.CATEGORIES[r.category]
        ax.axvline(r.x, color=colr, lw=0.4, alpha=0.6, zorder=1)
        ax.plot([r.x, r.x, xl], [1, 1 + 0.03 * f, 1 + 0.08 * f], transform=tr, color=colr, lw=0.6,
                clip_on=False)
        ax.plot(xl, 1 + 0.11 * f, transform=tr, marker=mk, ms=3, mfc=colr, mec="white", mew=0.3, ls="none",
                clip_on=False)
        ax.text(xl, 1 + 0.145 * f, style.italic_gene(r.label), transform=tr, rotation=90, fontsize=style.FS_TICK,
                ha="center", va="bottom", color=style.INK)


def main():
    style.setup()
    kg = key_genes_v4()
    print("Hu2021 VCF -> per-group allele counts")
    hu = hu_counts()
    hw = window_stats(hu, list(HU_GROUPS), HU_CONTRASTS, WIN_HU)
    hw["log2_He_sweet_over_field"] = np.log2(hw.he_sweet / hw.he_temperate_field)
    hw["log2_He_popcorn_over_field"] = np.log2(hw.he_popcorn / hw.he_temperate_field)
    hw.to_csv(TAB / "popgen_windows_Hu2021_1Mb.csv", index=False)

    print("Qu2026 array -> per-group allele counts")
    qu = qu_counts_shared_ref()
    for a, b in QU_CONTRASTS:
        ok = (qu[f"n_{a}"] >= MIN_ALLELES) & (qu[f"n_{b}"] >= MIN_ALLELES)
        pa, pb = qu[f"alt_{a}"] / qu[f"n_{a}"], qu[f"alt_{b}"] / qu[f"n_{b}"]
        num, den = hudson(pa, qu[f"n_{a}"], pb, qu[f"n_{b}"])
        qu[f"fst_{a}_vs_{b}"] = np.where(ok & (den > 0), num / den, np.nan)
    qu.to_csv(TAB / "popgen_snps_Qu2026.csv", index=False)
    qw = window_stats(qu, ["Sweet", "Waxy", "Field"], QU_CONTRASTS, WIN_QU)
    qw.to_csv(TAB / "popgen_windows_Qu2026_5Mb.csv", index=False)

    # statistics at each key gene: the window containing the gene start, and its genome-wide percentile
    recs = []
    for r in kg.dropna(subset=["chr"]).itertuples():
        rec = dict(symbol=r.symbol, category=r.category, v5_id=r.v5_id, v4_id=r.v4_id, chr=int(r.chr),
                   v4_start=int(r.v4_start))
        cols = ["fst_sweet_vs_temperate_field", "fst_popcorn_vs_temperate_field",
                "fst_sweet_vs_tropical", "log2_He_sweet_over_field"]
        row = hw[(hw.chr == r.chr) & (hw.win == int(r.v4_start) // WIN_HU)]
        for c in cols:   # Hu 2021: the 1-Mb window holding the gene, percentile among windows
            v = float(row[c].iloc[0]) if len(row) else np.nan
            rec[f"Hu_{c}"] = round(v, 4) if pd.notna(v) else np.nan
            rec[f"Hu_{c}_pct"] = round((hw[c] < v).mean() * 100, 1) if pd.notna(v) else np.nan
        near = qu[(qu.chr == r.chr) & ((qu.pos - r.v4_start).abs() <= QU_FLANK)]
        rec["Qu_n_snps_within_2Mb"] = len(near)
        for a, b in QU_CONTRASTS:  # Qu 2026: highest per-SNP Fst within +-2 Mb, percentile among SNPs
            c = f"fst_{a}_vs_{b}"
            v = near[c].max() if len(near) else np.nan
            rec[f"Qu_max_{c}"] = round(float(v), 4) if pd.notna(v) else np.nan
            rec[f"Qu_max_{c}_pct"] = round((qu[c] < v).mean() * 100, 1) if pd.notna(v) else np.nan
        recs.append(rec)
    kgs = pd.DataFrame(recs)
    kgs.to_csv(TAB / "popgen_key_genes.csv", index=False)
    print(kgs[["symbol", "Hu_fst_sweet_vs_temperate_field", "Hu_fst_sweet_vs_temperate_field_pct",
               "Hu_log2_He_sweet_over_field", "Qu_max_fst_Waxy_vs_Field", "Qu_max_fst_Waxy_vs_Field_pct",
               "Qu_max_fst_Sweet_vs_Field", "Qu_max_fst_Sweet_vs_Field_pct"]].to_string(index=False))

    # top 1% sweet-vs-field windows with the annotated genes inside
    genes4 = v4_gene_table()
    col = "fst_sweet_vs_temperate_field"
    thr = np.nanpercentile(hw[col], 99)
    top = hw[hw[col] >= thr].sort_values(col, ascending=False)
    rows = []
    for t in top.itertuples():
        inside = genes4[(genes4.chr == t.chr) & (genes4.end >= t.start) & (genes4.start < t.end)]
        named = inside[inside.locus_symbol.fillna("").ne("")]
        rows.append(dict(chr=t.chr, start_Mb=t.start / 1e6, end_Mb=t.end / 1e6, fst=round(getattr(t, col), 4),
                         log2_He_sweet_over_field=round(t.log2_He_sweet_over_field, 3),
                         diversity_reduced_in="sweet corn" if t.log2_He_sweet_over_field < 0 else "temperate field corn",
                         n_genes=len(inside),
                         named_genes=";".join(sorted(set(named.locus_symbol)))[:400],
                         key_genes=";".join(kg[(kg.chr == t.chr) & (kg.v4_start >= t.start) & (kg.v4_start < t.end)].symbol)))
    pd.DataFrame(rows).to_csv(TAB / "popgen_top_windows_sweet_vs_field.csv", index=False)

    # figure
    chrlen = hw.groupby("chr").end.max()
    offsets, acc = {}, 0
    for c in range(1, 11):
        offsets[c] = acc
        acc += chrlen[c] + 15e6
    W, H = style.WIDTH_IN, 7.9
    fig = plt.figure(figsize=(W, H))
    left, right, ah, zone, gap = 0.74, 0.36, 0.98, 0.38, 0.19   # inches; zone = gene labels above axes
    starch = ["sh2", "su1", "se1", "bt2", "sh1", "wx1", "ae1", "bt1", "du1", "su2"]
    fst = r"$F_\mathrm{ST}$"
    panels = [
        (hw, "fst_sweet_vs_temperate_field", fst + ", sweet vs\ntemperate field", starch, None),
        (hw, "log2_He_sweet_over_field", r"log$_\mathrm{2}$ $H_\mathrm{e}$ ratio," + "\nsweet / field", starch, 0),
        (hw, "fst_popcorn_vs_temperate_field", fst + ", popcorn vs\ntemperate field",
         ["p1", "ZmGa1", "tb1", "fl2", "o2"], None),
        (qu, "fst_Waxy_vs_Field", fst + ",\nwaxy vs field", ["wx1", "ae1", "su1", "sh2"], None),
        (qu, "fst_Sweet_vs_Field", fst + ",\nsweet vs field", ["sh2", "su1", "bt2", "se1", "wx1"], None),
    ]
    axes = []
    for k, (data, col, ylabel, genes, ref) in enumerate(panels):
        y0 = H - 0.03 - zone - ah - k * (zone + ah + gap)
        ax = style.axes_in(fig, left, y0, W - left - right, ah)
        manhattan(ax, data, col, kg, ylabel, genes, offsets, ref_line=ref)
        ax.set_xticks([(offsets[c] + chrlen[c] / 2) / 1e6 for c in range(1, 11)])
        ax.set_xticklabels([str(c) for c in range(1, 11)])
        style.panel_label(ax, "abcde"[k], dx=-(left - 0.04), dy=zone)
        axes.append(ax)
    axes[-1].set_xlabel("Chromosome (B73 RefGen_v4 coordinates)")
    style.save_png(fig, FIG / "Fig5_genome_scans_food_corn_selection")


if __name__ == "__main__":
    main()
