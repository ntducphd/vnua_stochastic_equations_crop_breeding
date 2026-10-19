"""
13_evidence_map.py -- evidence map of the records retained at screening: equation families and application domains.

Input: the records with final_decision = Include in pipeline_output/03_final_decisions.csv, joined to the Scopus and
Web of Science exports in raw/ for their author keywords and Scopus Index Keywords. The join repeats the matching
rule of step 01 (same normalised DOI, or same title reduced to letters and digits plus the same year), so every
retained record is linked to all the raw records it was merged from.

Each record is classified from its title, abstract, author keywords and Scopus Index Keywords with the regular
expressions in FAMILIES and DOMAINS. Both classifications are multi-label: a record carries every family and every
domain with at least one matching term, and a record without a domain term counts as general theory. Matching
ignores case and treats a space, a hyphen and every Unicode dash between words alike, so "Wright-Fisher",
"Wright–Fisher" and "Wright Fisher" are the same term.

Choices the term lists depend on, and why:
  - Web of Science Keywords Plus are not used. Clarivate generates them from the titles of the cited references, so
    they describe the reference list rather than the record (a theory paper citing Drosophila studies receives
    "DROSOPHILA-MELANOGASTER").
  - Scopus indexing terms that are assigned to whole classes of records are ignored: the check tags Human(s) and
    Animal(s), and the headings Markov chain(s), Stochastic control systems and Drug therapy (GENERIC_INDEX_TERMS).
  - "Markov chain" followed by "Monte Carlo" is an inference method, not a model of the system, and does not count.
  - "establishment" counts only as establishment probability or establishment of mutations, mutants, alleles,
    variants, clones or genes, or as a stand-alone keyword; the bare word also describes the establishment of
    clines, of viral latency or of plant colonies.
  - Games count for F5 only as mean-field or differential games. Evolutionary game dynamics are frequency-dependent
    selection in Moran or Wright-Fisher populations, with no controller or decision maker.
  - "individual-based" counts for F6 only when followed by "simulation"; an individual-based model that is analysed
    through its master equation or diffusion limit belongs to F1 or F2.
  - "Wright-Fisher controversy" and "debate" (the Fisher-Wright dispute), "Kingman" alone (Kingman's house-of-cards
    model and simplex), "Hamilton-Jacobi" alone (WKB and eikonal methods), "MDP" (maximal dissipation principle),
    "bird's eye" and the Scopus headings "Sexual Behavior, Animal", "Mating Preference, Animal" and "Animal
    Migration" do not count.
  - "introgression" is not a breeding term: in this set of records it describes natural hybridisation.
  - "seed", "leaf" and "leaves" are not plant terms ("lay the seeds for", "the leaves of the genealogy"); seed
    dormancy, dispersal and germination and seedlings are. "disease" preceded by "plant" is not a medical term.

Outputs:
  deposit/evidence_map_records.csv             one row per retained record: record_id (first Scopus EID, otherwise
                                               first Web of Science UT), all Scopus EIDs and WoS UTs, DOI, year,
                                               period, family flags F1-F6, no_family, domain flags; no titles or
                                               abstracts
  deposit/evidence_map_counts.csv              aggregate counts in long form (table, row, column, records):
                                               family x domain, family x period, plant-and-breeding x family
  ../figures/R01_fig_evidence_map.pdf, .png    a: records by family and domain; b: records per family by period

Usage: python -I pipeline/13_evidence_map.py   (run from data/, after step 03; the raw exports are
read as data only)
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True   # no __pycache__ folder in the deposited scripts folder
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, "scripts"))
from fig_style import apply, new_figure, panel, save, CYCLE, INK, GREY, MM  # noqa: E402
import matplotlib.colors as mcolors  # noqa: E402

FINAL_PATH = "pipeline_output/03_final_decisions.csv"
RECORDS_OUT = "deposit/evidence_map_records.csv"
COUNTS_OUT = "deposit/evidence_map_counts.csv"
FIG_DIR, FIG_STEM = "../figures", "R01_fig_evidence_map"

# ---------------------------------------------------------------------------------------------------------------
# Term lists
# ---------------------------------------------------------------------------------------------------------------
SEP = r"[\s\-‐‑‒–—―−]*"    # optional space, hyphen or Unicode dash
W2 = r"(?:[\w‐‑–—-]+{_}){0,2}"             # up to two intervening words


def terms(*parts):
    """One case-insensitive pattern from alternative terms; {_} stands for SEP."""
    return re.compile("|".join(p.replace("{_}", SEP) for p in parts), re.I)


FAMILY_NAMES = {
    "F1": "Master equations, Markov chains",
    "F2": "Diffusion approximations, SDEs",
    "F3": "Coalescent, genealogies",
    "F4": "Branching processes",
    "F5": "Stochastic control, decision",
    "F6": "Stochastic simulation",
}
FAMILIES = {
    "F1": terms(r"master{_}equations?", r"gillespie", r"stochastic{_}simulation{_}algorithms?", r"kinetic{_}monte{_}carlo",
                r"birth{_}(and{_})?death", r"\bmoran\b", r"markov{_}chains?(?!{_}monte)",
                r"(gene{_})?expression{_}noise", r"noise{_}in{_}(gene{_})?expression", r"stochastic{_}gene{_}expression",
                r"transcriptional{_}burst\w*", r"linear{_}noise{_}approximation", r"system{_}size{_}expansion"),
    "F2": terms(r"wright{_}fisher(?!{_}(controvers|debate))", r"fisher{_}wright(?!{_}(controvers|debate))",
                r"(?<!reaction)(?<!reaction-)(?<!reaction )(?<!reaction–)"
                r"diffusion{_}(approximations?|limits?|process(es)?|theory|equations?|models?)",
                r"jump{_}diffusions?", r"fokker{_}planck", r"kolmogorov(?!{_}smirnov)", r"kramers?{_}moyal",
                r"stochastic{_}(partial{_})?differential{_}equations?", r"ornstein{_}uhlenbeck", r"fixation{_}probabilit\w*",
                r"probabilit(y|ies){_}of{_}fixation", r"langevin", r"stratonovich", r"brownian{_}motion",
                r"wiener{_}process"),
    "F3": terms(r"coalescen\w*", r"genealog\w*", r"ancestral{_}(recombination|selection){_}graphs?", r"gene{_}trees?",
                r"most{_}recent{_}common{_}ancestors?"),
    "F4": terms(r"branching{_}process\w*", r"galton{_}watson", r"establishment{_}probabilit\w*",
                r"probabilit(y|ies){_}of{_}" + W2 + r"establishment",
                r"establishment{_}of{_}" + W2 + r"(mutations?|mutants?|alleles?|variants?|clones?|genes?)\b",
                r"(^|;)\s*establishment\s*(;|$)"),
    "F5": terms(r"optimal{_}control", r"stochastic{_}control", r"dynamic{_}programming", r"bellman",
                r"markov{_}decision", r"reinforcement{_}learning", r"mean{_}field{_}games?", r"differential{_}games?"),
    "F6": terms(r"individual{_}based{_}" + W2 + r"simulat\w*", r"agent{_}based", r"forward{_}" + W2 + r"simulat\w*",
                r"forward{_}time", r"simulation{_}stud(y|ies)", r"(?<!chain )monte{_}carlo{_}simulat\w*",
                r"exact{_}simulation", r"\bsimulators?\b",
                r"simulation{_}(programs?|software|tools?|packages?|frameworks?|platforms?)",
                r"\b(alphasimr?|qmsim|slim|simupop|mobps|xsim|qu{_}gene|quline|pedigreesim|fwdpp|quantinemo|msprime"
                r"|admixsim)\b"),
}
# Abbreviations are matched case-sensitively, in titles, abstracts and author keywords only.
FAMILY_ABBREVIATIONS = {"F2": re.compile(r"\bS(P)?DEs?\b"), "F5": re.compile(r"\bHJB\b")}

DOMAIN_NAMES = {
    "plant": "Plants, crops",
    "breeding": "Breeding",
    "animal": "Animals, livestock",
    "human": "Humans, medicine",
    "microbe": "Microbes, viruses, cancer",
    "none": "No domain term",
}
DOMAINS = {
    "plant": terms(r"\bplants?\b", r"\bcrops?\b", r"\bmaize\b", r"\bcorn\b", r"zea{_}mays", r"teosinte", r"\bwheat\b",
                   r"triticum", r"\brice\b", r"oryza", r"arabidopsis", r"\bcultivars?\b", r"\bbarley\b", r"sorghum",
                   r"soybeans?", r"\bcotton\b", r"tomato(es)?", r"potato(es)?", r"\bpeppers?\b", r"\bapples?\b", r"cassava",
                   r"brassica", r"rapeseed", r"canola", r"sunflower", r"millet", r"\boats?\b", r"\brye\b", r"grapevine",
                   r"tobacco", r"nicotiana", r"medicago", r"alfalfa", r"antirrhinum", r"snapdragon", r"mimulus",
                   r"eichhornia", r"seed{_}(dormancy|dispersal|germination)", r"seedlings?", r"pollen", r"pollinat\w*",
                   r"\bflower\w*", r"angiosperm\w*", r"gymnosperm\w*", r"conifer\w*", r"\bpines?\b", r"eucalypt\w*",
                   r"poplar", r"populus", r"\boaks?\b", r"quercus", r"agricultur\w*", r"agronom\w*", r"horticultur\w*",
                   r"germination", r"self{_}incompatib\w*"),
    "breeding": terms(r"(?<!timing of )(?<!skip )(?<!cooperative )\bbreed(ing|s|ers?)?\b"
                      r"(?!{_}(seasons?|sites?|grounds?|pairs?|success|systems?|habitats?))",
                      r"artificial{_}selection", r"genomic{_}selection", r"genomic{_}prediction", r"marker{_}assisted",
                      r"backcross\w*", r"genetic{_}gains?"),
    "animal": terms(r"\banimals?\b", r"livestock", r"cattle", r"\bcows?\b", r"dairy", r"\bbeef\b", r"\bpigs?\b", r"swine",
                    r"porcine", r"poultry", r"chickens?", r"\bsheep\b", r"\bgoats?\b", r"\bhorses?\b", r"\bdogs?\b",
                    r"\bbovine\b", r"\boxen\b", r"\bherds?\b", r"\bfish(es)?\b", r"fisher(y|ies)\b",
                    r"salmon(s|ids?|idae)?\b", r"trout", r"aquacult\w*", r"drosophila", r"\bflies\b",
                    r"(?<!on the )\bfly\b", r"\binsects?\b", r"mosquito(es)?", r"anopheles", r"aedes", r"whitefl(y|ies)",
                    r"aphids?", r"\bmoths?\b", r"tribolium", r"\bbees?\b", r"honeybees?", r"\bants?\b", r"butterfl\w*",
                    r"beetles?", r"\bmice\b", r"\bmouse\b", r"murine", r"\brats?\b", r"\bbirds?\b(?!['’]s{_}eye)", r"avian", r"mammal\w*",
                    r"primates?", r"\bvertebrates?\b", r"invertebrates?", r"nematodes?", r"elegans", r"zebrafish",
                    r"amphibians?", r"\bfrogs?\b", r"lizards?", r"snails?", r"daphnia", r"copepods?", r"\bwhales?\b",
                    r"\bdeer\b", r"\bwolf\b", r"\bwolves\b"),
    "human": terms(r"\bhumans?\b", r"patients?\b", r"clinic(al)?\b", r"\bmedic(al|ine|inal)\b", r"therap(y|ies|eutic)\w*",
                   r"hospital\w*", r"epidemiolog\w*", r"public{_}health", r"homo{_}sapiens", r"hominins?",
                   r"neanderthals?", r"(?<!plant )\bdiseases?\b"),
    "microbe": terms(r"bacteri\w*", r"microb\w*", r"\bvirus\w*", r"\bviral\b", r"\bvirions?\b", r"\b(bacterio)?phages?\b",
                     r"\byeasts?\b", r"e\.{_}coli", r"escherichia", r"saccharomyces", r"\bHIV\b", r"influenza", r"\bSARS\b",
                     r"pathogen\w*", r"plasmids?", r"\bfung(i|us|al)\b", r"archae(a|al|on)\b", r"micro{_}organisms?",
                     r"plasmodium", r"malaria", r"parasit\w*", r"cancer\w*", r"tumou?r\w*", r"carcinoma\w*",
                     r"oncogen\w*", r"leuka?emia", r"\bneoplas\w*", r"metasta(sis|ses|tic|si[sz]e\w*)\b"),
}
GENERIC_INDEX_TERMS = {"human", "humans", "animal", "animals", "markov chain", "markov chains",
                       "stochastic control systems", "drug therapy",
                       "sexual behavior, animal", "mating preference, animal", "animal migration"}

PERIODS = [("before 1980", 0, 1979), ("1980-1999", 1980, 1999), ("2000-2009", 2000, 2009),
           ("2010-2019", 2010, 2019), ("2020-2026", 2020, 2026)]
FAM = list(FAMILIES)
DOM = list(DOMAINS) + ["none"]


# ---------------------------------------------------------------------------------------------------------------
# Records and their keywords
# ---------------------------------------------------------------------------------------------------------------
def text(x):
    return "" if pd.isna(x) else str(x)


def norm_title(t):
    return re.sub(r"[^a-z0-9]", "", text(t).lower())


def norm_doi(d):
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", text(d).strip().lower())


def year4(y):
    m = re.search(r"(\d{4})", text(y))
    return m.group(1) if m else ""


def load_records():
    scopus = pd.read_csv(glob.glob("raw/scopus/*.csv")[0], encoding="utf-8-sig")
    wos = pd.read_excel(glob.glob("raw/WoS/*.xls")[0])
    raw = pd.concat([
        pd.DataFrame({"db": "Scopus", "rid": scopus["EID"], "title": scopus["Title"], "year": scopus["Year"],
                      "doi": scopus["DOI"], "kw_author": scopus["Author Keywords"],
                      "kw_index": scopus["Index Keywords"]}),
        pd.DataFrame({"db": "Web of Science", "rid": wos["UT (Unique WOS ID)"], "title": wos["Article Title"],
                      "year": wos["Publication Year"], "doi": wos["DOI"], "kw_author": wos["Author Keywords"],
                      "kw_index": np.nan}),
    ], ignore_index=True)
    raw["tk"] = raw["title"].map(norm_title) + "|" + raw["year"].map(year4)
    raw["dk"] = raw["doi"].map(norm_doi)

    # union-find over the two signals of step 01
    parent = list(range(len(raw)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for key, usable in (("dk", raw["dk"] != ""), ("tk", raw["title"].map(norm_title) != "")):
        groups = {}
        for i in raw.index[usable]:
            groups.setdefault(raw.at[i, key], []).append(i)
        for members in groups.values():
            for j in members[1:]:
                ri, rj = find(members[0]), find(j)
                if ri != rj:
                    parent[rj] = ri
    raw["group"] = [find(i) for i in raw.index]

    final = pd.read_csv(FINAL_PATH, encoding="utf-8-sig")
    final["group"] = (final["title"].map(norm_title) + "|" + final["year"].map(year4)).map(
        dict(zip(raw["tk"], raw["group"])))
    assert final["group"].notna().all() and final["group"].is_unique, "screened records do not match the raw exports"
    assert set(final["group"]) == set(raw["group"]), "raw records left without a screened record"

    def joined(series):
        return "; ".join(text(x) for x in series if text(x))

    agg = raw.groupby("group").apply(lambda g: pd.Series({
        "scopus_eid": ";".join(sorted(g.loc[g["db"] == "Scopus", "rid"].astype(str))),
        "wos_ut": ";".join(sorted(g.loc[g["db"] == "Web of Science", "rid"].astype(str))),
        "kw_author": joined(g["kw_author"]),
        "kw_index": joined(g["kw_index"]),
    }), include_groups=False)
    final = final.merge(agg, left_on="group", right_index=True, how="left")
    retained = final[final["final_decision"] == "Include"].reset_index(drop=True)
    retained["record_id"] = [e.split(";")[0] if e else u.split(";")[0]
                             for e, u in zip(retained["scopus_eid"], retained["wos_ut"])]
    assert retained["record_id"].is_unique
    return retained


# ---------------------------------------------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------------------------------------------
def classify(rec):
    abstract = text(rec["abstract"])
    if abstract.strip().lower() == "[no abstract available]":
        abstract = ""
    ta = text(rec["title"]) + " | " + abstract
    ka = text(rec["kw_author"])
    ki = "; ".join(k for k in (x.strip() for x in text(rec["kw_index"]).split(";"))
                   if k and k.lower() not in GENERIC_INDEX_TERMS)
    out = {}
    for code, pattern in FAMILIES.items():
        hit = any(pattern.search(x) for x in (ta, ka, ki))
        if code in FAMILY_ABBREVIATIONS:
            hit = hit or any(FAMILY_ABBREVIATIONS[code].search(x) for x in (ta, ka))
        out[code] = int(hit)
    out["no_family"] = int(not any(out[c] for c in FAM))
    for code, pattern in DOMAINS.items():
        out[code] = int(any(pattern.search(x) for x in (ta, ka, ki)))
    out["none"] = int(not any(out[c] for c in DOMAINS))
    return out


def period(y):
    y = int(year4(y))
    return next(label for label, lo, hi in PERIODS if lo <= y <= hi)


# ---------------------------------------------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------------------------------------------
def draw_figure(fd, fp, n_family, n_domain, n_period):
    apply()
    W, H = 160.0, 112.0
    fig, _ax = new_figure(height_mm=H, constrained_layout=False)
    _ax.remove()

    def axes_mm(left, bottom, width, height):
        return fig.add_axes([left / W, bottom / H, width / W, height / H])

    # a: heat map, families (rows) x domains (columns); shading follows the square root of the count
    rows = FAM + ["no_family"]
    cw, rh = 14.0, 6.0
    left_a, bottom_a = 44.0, 56.0
    ax = axes_mm(left_a, bottom_a, cw * len(DOM), rh * len(rows))
    counts = fd.loc[rows, DOM].values.astype(float)
    cmap = mcolors.LinearSegmentedColormap.from_list("oi_blue", ["#FFFFFF", CYCLE[0]])
    norm = mcolors.PowerNorm(gamma=0.5, vmin=0, vmax=counts.max())
    ax.pcolormesh(np.arange(len(DOM) + 1), np.arange(len(rows) + 1), counts, cmap=cmap, norm=norm,
                  edgecolors="white", linewidth=0.8)
    for i in range(len(rows)):
        for j in range(len(DOM)):
            v = counts[i, j]
            ax.text(j + 0.5, i + 0.5, f"{int(v)}", ha="center", va="center",
                    color="white" if norm(v) > 0.6 else (GREY if v == 0 else INK))
    ax.axhline(len(FAM), color="white", lw=2.5)
    ax.axvline(len(DOM) - 1, color="white", lw=2.5)
    ax.set_xlim(0, len(DOM))
    ax.set_ylim(len(rows), 0)
    ax.set_yticks(np.arange(len(rows)) + 0.5)
    ax.set_yticklabels([FAMILY_NAMES.get(r, "No family term") for r in rows])
    ax.set_xticks(np.arange(len(DOM)) + 0.5)
    ax.set_xticklabels([DOMAIN_NAMES[d].replace(", ", ",\n").replace("No domain term", "No domain\nterm")
                        + f"\n({n_domain[d]})" for d in DOM], linespacing=1.1)
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=2)
    for sp in ax.spines.values():
        sp.set_visible(False)
    x_all = len(DOM) + 0.5
    ax.text(x_all, -0.15, "All", ha="center", va="bottom")
    for i, r in enumerate(rows):
        ax.text(x_all, i + 0.5, f"{n_family[r]}", ha="center", va="center")
    panel(ax, "a", x=(3.5 - left_a) / (cw * len(DOM)), y=(H - 4.0 - bottom_a) / (rh * len(rows)))

    # b: records per family by publication period
    left_b, bottom_b, width_b, height_b = 17.0, 11.0, 90.0, 36.0
    bx = axes_mm(left_b, bottom_b, width_b, height_b)
    markers = ["o", "s", "^", "D", "v", "P"]
    x = np.arange(len(PERIODS))
    for k, c in enumerate(FAM):
        bx.plot(x, fp.loc[c, [p for p, _, _ in PERIODS]].values, color=CYCLE[k], marker=markers[k], markersize=3.2,
                label=FAMILY_NAMES[c], clip_on=False, zorder=3)
    bx.set_xticks(x)
    bx.set_xticklabels([p.replace("-", "\u2013") + f"\n({n_period[p]})" for p, _, _ in PERIODS], linespacing=1.1)
    bx.set_xlim(-0.3, len(PERIODS) - 0.7)
    bx.set_ylim(0, None)
    bx.set_ylabel("Records")
    bx.set_xlabel("Publication period (retained records)")
    bx.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False, handlelength=1.8)
    panel(bx, "b", x=(3.5 - left_b) / width_b, y=(51.0 - bottom_b) / height_b)
    return save(fig, FIG_DIR, FIG_STEM)


# ---------------------------------------------------------------------------------------------------------------
def main():
    os.makedirs("deposit", exist_ok=True)
    retained = load_records()
    flags = pd.DataFrame([classify(r) for _, r in retained.iterrows()])
    out = pd.concat([retained[["record_id", "scopus_eid", "wos_ut", "doi", "year"]], flags], axis=1)
    out["period"] = retained["year"].map(period)
    out = out[["record_id", "scopus_eid", "wos_ut", "doi", "year", "period"] + FAM + ["no_family"] + DOM]
    out.to_csv(RECORDS_OUT, index=False, encoding="utf-8-sig")

    rows = FAM + ["no_family"]
    fd = pd.DataFrame({d: [int(((out[r] == 1) & (out[d] == 1)).sum()) for r in rows] for d in DOM}, index=rows)
    fp = pd.DataFrame({p: [int(((out[r] == 1) & (out["period"] == p)).sum()) for r in rows] for p, _, _ in PERIODS},
                      index=rows)
    n_family = {r: int(out[r].sum()) for r in rows}
    n_domain = {d: int(out[d].sum()) for d in DOM}
    n_period = {p: int((out["period"] == p).sum()) for p, _, _ in PERIODS}
    plant_breeding = (out["plant"] == 1) & (out["breeding"] == 1)
    pb = {r: int((plant_breeding & (out[r] == 1)).sum()) for r in rows}

    long = []
    for r in rows:
        long += [("family_x_domain", r, d, fd.at[r, d]) for d in DOM]
        long.append(("family_x_domain", r, "all", n_family[r]))
    long += [("family_x_domain", "all", d, n_domain[d]) for d in DOM]
    long.append(("family_x_domain", "all", "all", len(out)))
    for r in rows:
        long += [("family_x_period", r, p, fp.at[r, p]) for p, _, _ in PERIODS]
    long += [("family_x_period", "all", p, n_period[p]) for p, _, _ in PERIODS]
    long += [("plant_and_breeding_x_family", r, "plant_and_breeding", pb[r]) for r in rows]
    long.append(("plant_and_breeding_x_family", "all", "plant_and_breeding", int(plant_breeding.sum())))
    pd.DataFrame(long, columns=["table", "row", "column", "records"]).to_csv(COUNTS_OUT, index=False,
                                                                               encoding="utf-8-sig")
    png = draw_figure(fd, fp, n_family, n_domain, n_period)

    pd.set_option("display.width", 160)
    print(f"retained records: {len(out)}")
    print("\nfamily x domain (multi-label; 'all' = records with the family):")
    print(fd.assign(all=pd.Series(n_family)).to_string())
    print("records per domain:", n_domain)
    print("\nfamily x publication period:")
    print(fp.to_string())
    print("records per period:", n_period)
    print(f"\nrecords with no family term: {n_family['no_family']}")
    print(f"records with a plant/crop term and a breeding term: {int(plant_breeding.sum())}; with each family: {pb}")
    print("record_id of those:", ", ".join(out.loc[plant_breeding, "record_id"]))
    print(f"\nSaved {RECORDS_OUT}, {COUNTS_OUT} and {png} (+ PDF)")


if __name__ == "__main__":
    main()
