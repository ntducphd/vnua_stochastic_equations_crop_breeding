"""
12_build_deposit_tables.py -- the redistributable record of the literature search, written to deposit/.

The Scopus and Web of Science exports in raw/ and the pipeline_output/ CSV files derived from them carry publisher
abstracts and author contact details and cannot be redistributed under the providers' terms. This step writes
what can be deposited instead:

  deposit/records_identified.csv       all 1,695 retrieved records: database, Scopus EID or Web of Science UT,
                                       DOI, title, year, source title, document type (no abstracts, no addresses)
  deposit/screening_decisions.csv      one row per screened (deduplicated) record with the classifier suggestion
                                       and the final decision, reason and notes; reviewer_source is
                                       ai_confirmed_default (the AI assistant kept the classifier suggestion) or
                                       ai_ledger_decision (decided through a ledger entry); no abstracts
  deposit/screening_decisions_ledger.csv  the 202 ledger decisions, each with its reason
  deposit/R01_screening_workbook.xlsx  the workbook of step 05 (no abstracts)
  deposit/citation_reconciliation.csv  the works cited in the manuscript matched against the screened records
                                       (step 04)

All screening decisions were made by the AI assistant (Claude, Anthropic) at the first author's instruction; the
first author takes responsibility for them (manuscript Section 2.2). The search terms are given in manuscript
Section 2.1; the exact query strings were not retained. The record-level calibration readings behind kappa were
not kept; only the summary table survives (workbook sheet Calibration_Validation).

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
assert set(final["reviewer_source"]) <= {"ai_confirmed_default", "ai_ledger_decision"}
final.drop(columns=["abstract"]).to_csv(f"{OUT}/screening_decisions.csv", index=False, encoding="utf-8-sig")

# The deposited decisions must cover every screened record (the deduplicated output of step 01), each once.
KEY = ["title", "year", "doi", "database"]


def record_keys(d):
    return sorted("|".join("" if pd.isna(v) else str(v) for v in row) for row in d[KEY].itertuples(index=False))


deposited = pd.read_csv(f"{OUT}/screening_decisions.csv", encoding="utf-8-sig")
assert record_keys(deposited) == record_keys(screened), "screening_decisions.csv does not cover every screened record"
assert deposited["final_decision"].isin(["Include", "Maybe", "Exclude"]).all()

shutil.copyfile("pipeline/screening_decisions_ledger.csv", f"{OUT}/screening_decisions_ledger.csv")
shutil.copyfile("pipeline_output/R01_PRISMA_screening_FINAL.xlsx", f"{OUT}/R01_screening_workbook.xlsx")
shutil.copyfile("pipeline_output/04_citation_reconciliation.csv", f"{OUT}/citation_reconciliation.csv")

changed = int((final["final_decision"] != final["ai_suggested_decision"]).sum())
recon = pd.read_csv(f"{OUT}/citation_reconciliation.csv", encoding="utf-8-sig")
print("records identified: %d (Scopus %d, Web of Science %d)" % (len(ids), len(scopus), len(wos)))
print("screened: %d; final %s; changed by the AI assistant: %d"
      % (len(final), final["final_decision"].value_counts().to_dict(), changed))
print("cited works: %d, also retrieved by the search: %d" % (len(recon), int(recon["found_in_search"].sum())))
print("written to %s/" % OUT)
