"""
03_merge_screening_decisions.py -- final screening decisions: the classifier suggestions of step 02, confirmed or
overridden through the decision ledger.

At the first author's instruction, the AI assistant (Claude, Anthropic) confirmed or overrode every classifier
suggestion; the first author takes responsibility for the decisions (manuscript Section 2.2):
  - all 15 classifier-Exclude records (9 changed: 7 to Include, 2 to Maybe)
  - all 157 classifier-Maybe records (each decided individually)
  - the medium-confidence classifier-Include records, record by record
  - the high-confidence classifier-Include records, checked from their titles in batches
pipeline/screening_decisions_ledger.csv records each decision made through a ledger entry, with its reason: every
decision that differs from the classifier suggestion, every decision on a classifier-Maybe record, and the
confirmation of six classifier-Exclude suggestions. A record with no ledger entry keeps the classifier suggestion
(reviewer_source = ai_confirmed_default).

Ledger rows are matched to records by title. Both titles are reduced to their letters and digits after Unicode
compatibility normalisation and case folding, so every dash (hyphen, en dash, em dash, minus sign), apostrophe,
quotation mark and other punctuation mark, and all spacing, is treated the same way. The script stops with an
error when a ledger row matches no record or more than one, so that no ledger decision is dropped or applied to
the wrong record.

Usage: python -I pipeline/03_merge_screening_decisions.py   (run from data/)
"""
import sys
import unicodedata

import pandas as pd

SCREENED_PATH = "pipeline_output/02_merged_with_ai_screening.csv"
LEDGER_PATH = "pipeline/screening_decisions_ledger.csv"
OUT_PATH = "pipeline_output/03_final_decisions.csv"


def norm_title(t):
    """Matching key: NFKC-normalised, case-folded title reduced to its letters and digits."""
    if not isinstance(t, str):
        return ""
    t = unicodedata.normalize("NFKC", t).casefold()
    return "".join(ch for ch in t if ch.isalnum())


screened = pd.read_csv(SCREENED_PATH, encoding="utf-8-sig")
screened["_norm_title"] = screened["title"].apply(norm_title)

ledger = pd.read_csv(LEDGER_PATH, encoding="utf-8-sig")
ledger["_norm_title"] = ledger["title"].apply(norm_title)

# Every ledger row must match exactly one screened record.
n_records = screened["_norm_title"].value_counts()
ledger["_n_matches"] = ledger["_norm_title"].map(n_records).fillna(0).astype(int)
unmatched = ledger[ledger["_n_matches"] != 1]
if len(unmatched):
    for _, e in unmatched.iterrows():
        print(f"ERROR: ledger row matches {e['_n_matches']} screened records: {e['title']}")
    sys.exit(f"{len(unmatched)} ledger row(s) do not match exactly one screened record; no output written.")

dupes = ledger[ledger["_norm_title"].duplicated(keep=False)].sort_values("_norm_title")
if len(dupes):
    print("WARNING duplicate titles in ledger:")
    print(dupes)

ledger_map = ledger.drop_duplicates(subset="_norm_title", keep="last").set_index("_norm_title")


def resolve(row):
    t = row["_norm_title"]
    if t in ledger_map.index:
        e = ledger_map.loc[t]
        return pd.Series({
            "final_decision": e["final_decision"],
            "final_exclusion_reason": e["final_exclusion_reason"] if pd.notna(e["final_exclusion_reason"]) else "",
            "final_notes": e["final_notes"] if pd.notna(e["final_notes"]) else "",
            "reviewer_source": "ai_ledger_decision",
        })
    ai = row["ai_suggested_decision"]
    if ai == "Maybe":
        return pd.Series({
            "final_decision": "NEEDS_REVIEW",
            "final_exclusion_reason": "",
            "final_notes": "no ledger decision found for a classifier-Maybe record - gap, must be resolved before finalising",
            "reviewer_source": "GAP",
        })
    return pd.Series({
        "final_decision": ai,
        "final_exclusion_reason": "",
        "final_notes": "confirmed by the AI assistant - matches the classifier suggestion, no exception raised",
        "reviewer_source": "ai_confirmed_default",
    })


resolved = screened.apply(resolve, axis=1)
final = pd.concat([screened, resolved], axis=1).drop(columns=["_norm_title"])

print("\n--- final_decision counts ---")
print(final["final_decision"].value_counts())
print("\n--- reviewer_source counts ---")
print(final["reviewer_source"].value_counts())

gaps = final[final["final_decision"] == "NEEDS_REVIEW"]
if len(gaps):
    print(f"\nERROR: {len(gaps)} unresolved AI-Maybe records with no ledger entry:")
    print(gaps[["title", "ai_suggested_decision"]].to_string())
else:
    print("\nNo gaps: every record has an explicit final decision.")

final.to_csv(OUT_PATH, index=False, encoding="utf-8-sig")
print(f"\nSaved {OUT_PATH}")
