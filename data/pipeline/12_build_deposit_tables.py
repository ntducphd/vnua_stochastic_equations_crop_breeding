"""
12_build_deposit_tables.py -- the redistributable record of the literature search, written to deposit/.

The Scopus and Web of Science exports in raw/ and the pipeline_output/ CSV files derived from them carry publisher
abstracts, keywords and author contact details and cannot be redistributed under the providers' terms. This step
writes what can be deposited instead:

  deposit/records_identified.csv       all 1,695 retrieved records: database, Scopus EID or Web of Science UT,
                                       DOI, title, year, source title, document type (no abstracts, no addresses)
  deposit/screening_decisions.csv      one row per screened (deduplicated) record with the decision, the reason,
                                       the rule that decided it, its scores and matched terms, and the version and
                                       SHA-256 of the rules file (no abstracts, no keywords)
  deposit/R01_screening_workbook.xlsx  the workbook of step 05
  deposit/citation_reconciliation.csv  the works cited in the manuscript matched against the screened records
                                       (step 04)
  deposit/screening_check.csv          written by pipeline/03_check_screening.py score: the blinded reading of the
                                       check sample next to the rules' decision

The decisions are those of the screening rules (pipeline/02_screen_records.py, pipeline/screening_rules.yaml;
manuscript Section 2.2). The search terms are given in manuscript Section 2.1; the exact query strings were not
retained.

Usage: python -I pipeline/12_build_deposit_tables.py   (run from data/, after steps 01-05;
the raw exports are read as data only)
"""
import glob
import os
import shutil

import pandas as pd

OUT = "deposit"
os.makedirs(OUT, exist_ok=True)

scopus = pd.read_csv(glob.glob("raw/scopus/*.csv")[0], encoding="utf-8-sig")
wos = pd.read_excel(glob.glob("raw/WoS/*.xls")[0])
ids = pd.concat([
    pd.DataFrame({"database": "Scopus", "record_id": scopus["EID"], "doi": scopus["DOI"], "title": scopus["Title"],
                  "year": scopus["Year"], "source_title": scopus["Source title"],
                  "document_type": scopus["Document Type"]}),
    pd.DataFrame({"database": "Web of Science", "record_id": wos["UT (Unique WOS ID)"], "doi": wos["DOI"],
                  "title": wos["Article Title"], "year": wos["Publication Year"], "source_title": wos["Source Title"],
                  "document_type": wos["Document Type"]}),
], ignore_index=True)
assert len(ids) == 1695, len(ids)
ids.to_csv(f"{OUT}/records_identified.csv", index=False, encoding="utf-8-sig")

screened = pd.read_csv("pipeline_output/01_merged_deduplicated.csv", encoding="utf-8-sig")
final = pd.read_csv("pipeline_output/03_final_decisions.csv", encoding="utf-8-sig")
assert set(final["reviewer_source"]) == {"screening_rules"}
final.drop(columns=["abstract", "keywords", "context_terms"]).to_csv(
    f"{OUT}/screening_decisions.csv", index=False, encoding="utf-8-sig")

# The deposited decisions must cover every screened record (the deduplicated output of step 01), each once.
KEY = ["title", "year", "doi", "database"]


def record_keys(d):
    return sorted("|".join("" if pd.isna(v) else str(v) for v in row) for row in d[KEY].itertuples(index=False))


deposited = pd.read_csv(f"{OUT}/screening_decisions.csv", encoding="utf-8-sig")
assert record_keys(deposited) == record_keys(screened), "screening_decisions.csv does not cover every screened record"
assert deposited["final_decision"].isin(["Include", "Maybe", "Exclude"]).all()

shutil.copyfile("pipeline_output/R01_PRISMA_screening_FINAL.xlsx", f"{OUT}/R01_screening_workbook.xlsx")
shutil.copyfile("pipeline_output/04_citation_reconciliation.csv", f"{OUT}/citation_reconciliation.csv")

recon = pd.read_csv(f"{OUT}/citation_reconciliation.csv", encoding="utf-8-sig")
print("records identified: %d (Scopus %d, Web of Science %d)" % (len(ids), len(scopus), len(wos)))
print("screened: %d; decisions %s; rules %s" % (len(final), final["final_decision"].value_counts().to_dict(),
                                               final["rules_version"].iloc[0]))
print("cited works: %d, also retrieved by the search: %d" % (len(recon), int(recon["found_in_search"].sum())))
print("written to %s/" % OUT)
