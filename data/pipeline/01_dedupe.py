"""
01_dedupe.py -- PRISMA identification stage: merge the Scopus and Web of Science exports and remove duplicates.

Two records are taken to be the same work when their normalised DOIs match OR their normalised titles and
publication years match. The two signals are evaluated independently and joined with a union-find, so a record
with a DOI still merges with a copy that has none, and a work indexed under two different DOIs (older Applied
Probability Trust articles carry both a Project Euclid and a Cambridge Core DOI) still merges through its title
and year. Titles are compared in lower case with every character other than a-z and 0-9 removed, so copies that
differ only in spacing, hyphenation or punctuation ("OVER-DOMINANT" and "overdominant") give the same key.
The script reads only the raw exports in raw/.

Usage: python -I pipeline/01_dedupe.py   (run from data/)
"""
import re
import pandas as pd

SCOPUS_PATH = "raw/scopus/scopus_export_Jul 26-2026_787a0001-ab41-4be3-8cb8-054f18dc6721.csv"
WOS_PATH = "raw/WoS/savedrecs.xls"
OUT_PATH = "pipeline_output/01_merged_deduplicated.csv"
SUMMARY_PATH = "pipeline_output/01_prisma_identification_summary.txt"


def norm_title(t):
    """Title key: lower case, letters a-z and digits only (spaces, hyphens, dashes and punctuation removed)."""
    if pd.isna(t):
        return ""
    return re.sub(r"[^a-z0-9]", "", str(t).lower())


def norm_doi(d):
    if pd.isna(d):
        return ""
    d = str(d).strip().lower()
    d = re.sub(r"^https?://(dx\.)?doi\.org/", "", d)
    return d


# ---- load Scopus ----
scopus = pd.read_csv(SCOPUS_PATH, encoding="utf-8-sig")
scopus_std = pd.DataFrame({
    "title": scopus["Title"],
    "authors": scopus["Authors"],
    "year": scopus["Year"],
    "source_title": scopus["Source title"],
    "doi": scopus["DOI"],
    "abstract": scopus.get("Abstract", ""),
    "doc_type": scopus.get("Document Type", ""),
    "database": "Scopus",
})

# ---- load Web of Science ----
wos = pd.read_excel(WOS_PATH)
wos_std = pd.DataFrame({
    "title": wos["Article Title"],
    "authors": wos["Authors"],
    "year": wos["Publication Year"],
    "source_title": wos["Source Title"],
    "doi": wos["DOI"],
    "abstract": wos.get("Abstract", ""),
    "doc_type": wos.get("Document Type", ""),
    "database": "Web of Science",
})

n_scopus_raw = len(scopus_std)
n_wos_raw = len(wos_std)

combined = pd.concat([scopus_std, wos_std], ignore_index=True).reset_index(drop=True)
combined["_doi_norm"] = combined["doi"].apply(norm_doi)
combined["_title_norm"] = combined["title"].apply(norm_title)
combined["_year_str"] = combined["year"].astype(str).str.extract(r"(\d{4})")[0]
combined["_ty_key"] = combined["_title_norm"] + "|" + combined["_year_str"].astype(str)

n = len(combined)

# ---- union-find over BOTH signals independently (DOI match, title+year match) ----
parent = list(range(n))


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra


by_doi = {}
for i, d in enumerate(combined["_doi_norm"]):
    if d:
        by_doi.setdefault(d, []).append(i)
for idxs in by_doi.values():
    for j in idxs[1:]:
        union(idxs[0], j)

by_ty = {}
for i, k in enumerate(combined["_ty_key"]):
    if combined.loc[i, "_title_norm"]:  # guard against blank titles colliding
        by_ty.setdefault(k, []).append(i)
for idxs in by_ty.values():
    for j in idxs[1:]:
        union(idxs[0], j)

combined["_group"] = [find(i) for i in range(n)]

n_found_both_or_dup = combined["_group"].duplicated().sum()  # rows that are NOT the first of their group

dedup_rows = []
for gid, group in combined.groupby("_group"):
    dbs = sorted(set(group["database"]))
    # Prefer a representative row whose title is NOT all-uppercase (WoS legacy
    # exports render many older titles in ALL CAPS; Scopus usually keeps
    # normal title case for the same paper), then a DOI present, then the
    # Scopus copy of the abstract over the WoS copy (WoS's legacy .xls export
    # of this corpus was found to break hyphens/subscripts with stray spaces,
    # e.g. "Wright- Fisher" and "A 1" instead of "Wright-Fisher" and "A1" -
    # inflating length while corrupting the exact-substring term matching the
    # screening algorithm relies on), then the longest abstract as a final
    # tiebreak when neither source is Scopus.
    group = group.copy()
    title_str = group["title"].fillna("").astype(str)
    group["_not_allcaps"] = ~(title_str.str.upper() == title_str) | (title_str.str.len() == 0)
    group["_has_doi"] = group["_doi_norm"] != ""
    group["_is_scopus"] = group["database"] == "Scopus"
    group["_abs_len"] = group["abstract"].fillna("").astype(str).str.len()
    rep = group.sort_values(
        ["_not_allcaps", "_has_doi", "_is_scopus", "_abs_len"], ascending=[False, False, False, False]
    ).iloc[0].copy()
    rep["database"] = " + ".join(dbs)
    dedup_rows.append(rep)

deduped = pd.DataFrame(dedup_rows).drop(
    columns=["_doi_norm", "_title_norm", "_year_str", "_ty_key", "_group",
             "_not_allcaps", "_has_doi", "_is_scopus", "_abs_len"]
)
deduped = deduped.sort_values(["year", "title"], ascending=[False, True]).reset_index(drop=True)
deduped.to_csv(OUT_PATH, index=False, encoding="utf-8-sig")

n_after_dedup = len(deduped)
n_removed = (n_scopus_raw + n_wos_raw) - n_after_dedup
n_multi_source_groups = sum(1 for gid, group in combined.groupby("_group") if len(set(group["database"])) > 1)

summary = f"""PRISMA identification stage ({pd.Timestamp.now().strftime('%Y-%m-%d')})

Database searches:
  Scopus:          {n_scopus_raw} records
  Web of Science:  {n_wos_raw} records
  ---------------------------------
  Total identified: {n_scopus_raw + n_wos_raw} records

Deduplication (union-find over a DOI match OR a title + year match, evaluated
independently; titles compared in lower case on letters and digits only):
  Records found in more than one database: {n_multi_source_groups}
  Duplicate records removed: {n_removed}
  ---------------------------------
  Unique records after deduplication: {n_after_dedup}

Output: {OUT_PATH} ({n_after_dedup} rows, input to the title/abstract screening)
"""
with open(SUMMARY_PATH, "w", encoding="utf-8") as f:
    f.write(summary)

print(summary)
