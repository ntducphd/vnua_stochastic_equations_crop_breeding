"""
04_reconcile_citations.py -- cross-check the references cited in the manuscript (data/cited_references.csv, one
row per numbered entry of the reference list, exported from the manuscript of 13 October 2026) against the
screened records. A cited work counts as retrieved when its DOI, or failing that its normalised title, matches a
screened record. The search was not used to choose the cited works (manuscript Section 2.3).

Usage: python pipeline/04_reconcile_citations.py   (run from data/)
"""
import re
import pandas as pd

CITED_PATH = "cited_references.csv"
SCREENED_PATH = "pipeline_output/03_final_decisions.csv"
OUT_PATH = "pipeline_output/04_citation_reconciliation.csv"


def norm_doi(d):
    if pd.isna(d):
        return ""
    d = str(d).strip().lower()
    d = re.sub(r"^https?://(dx\.)?doi\.org/", "", d)
    return d


def norm_title(t):
    if pd.isna(t):
        return ""
    t = str(t).lower()
    t = re.sub(r"[^a-z0-9 ]", " ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


cited = pd.read_csv(CITED_PATH, encoding="utf-8-sig")
screened = pd.read_csv(SCREENED_PATH, encoding="utf-8-sig")

cited["_doi_norm"] = cited["DOI"].apply(norm_doi)
cited["_title_norm"] = cited["Reference"].apply(norm_title)
screened["_doi_norm"] = screened["doi"].apply(norm_doi)
screened["_title_norm"] = screened["title"].apply(norm_title)

doi_map = screened.set_index("_doi_norm")
title_map = screened.drop_duplicates(subset="_title_norm").set_index("_title_norm")

rows = []
for _, r in cited.iterrows():
    match = None
    match_method = None
    if r["_doi_norm"] and r["_doi_norm"] in doi_map.index:
        m = doi_map.loc[r["_doi_norm"]]
        if isinstance(m, pd.DataFrame):
            m = m.iloc[0]
        match = m
        match_method = "DOI"
    else:
        # title fallback: a screened title (at least 25 characters) contained in the reference entry
        hits = [t for t in title_map.index if len(t) >= 25 and t in r["_title_norm"]]
        if hits:
            match = title_map.loc[max(hits, key=len)]
            match_method = "title"

    if match is not None:
        rows.append({
            "ManuscriptRefNum": r["ManuscriptRefNum"], "FirstAuthor": r["FirstAuthor"], "Year": r["Year"],
            "Reference": r["Reference"], "found_in_search": True, "match_method": match_method,
            "screening_decision": match["final_decision"],
        })
    else:
        rows.append({
            "ManuscriptRefNum": r["ManuscriptRefNum"], "FirstAuthor": r["FirstAuthor"], "Year": r["Year"],
            "Reference": r["Reference"], "found_in_search": False, "match_method": "",
            "screening_decision": "not_found_by_search",
        })

out = pd.DataFrame(rows)
out.to_csv(OUT_PATH, index=False, encoding="utf-8-sig")

n_total = len(out)
n_found = out["found_in_search"].sum()
n_not_found = n_total - n_found
print(f"Currently-cited references: {n_total}")
print(f"  Found by the Scopus/WoS search: {n_found}")
print(f"  NOT found by the search (predate/outside the Boolean query): {n_not_found}")
print()
print("Of those found, screening decision breakdown:")
print(out[out["found_in_search"]]["screening_decision"].value_counts())

conflicts = out[(out["found_in_search"]) & (out["screening_decision"] != "Include")]
if len(conflicts):
    print(f"\n{len(conflicts)} currently-cited references were found by the search but NOT screened Include - real conflicts to review:")
    print(conflicts[["ManuscriptRefNum", "FirstAuthor", "Year", "screening_decision"]].to_string())

not_found = out[~out["found_in_search"]]
print(f"\n{len(not_found)} currently-cited references NOT found by the search:")
print(not_found[["ManuscriptRefNum", "FirstAuthor", "Year"]].to_string())

print(f"\nSaved {OUT_PATH}")
