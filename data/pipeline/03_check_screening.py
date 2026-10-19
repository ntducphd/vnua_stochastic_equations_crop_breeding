"""
03_check_screening.py -- check the screening rules against a blinded human reading of a random sample.

  python -I pipeline/03_check_screening.py sample   draw the check sample and write the reading sheet
  python -I pipeline/03_check_screening.py short    shorten the reading sheet to SHORT_SIZE per class (see below)
  python -I pipeline/03_check_screening.py score    score the filled sheet against the rules' decisions

sample: draws, with a fixed seed, a sample stratified by the rules' decision (SAMPLE_SIZE per class, or the
whole class when it is smaller) and writes
  check/screening_check_sheet.xlsx   title, abstract, keywords, source and document type of each sampled record,
                                     in random order, with empty columns for the reader's decision and reason;
                                     the rules' decisions are NOT in this file
  check/screening_check_key.csv      which record each check_id is and the rules' decision (do not open it
                                     before the sheet is filled)
The reader applies ../protocol/eligibility_criteria.md to the title, abstract and keywords only.

short: keeps, from the sample in the key, SHORT_SIZE records per class (seed SEED + 1), leaving out records marked
revealed_before_reading, and rewrites the sheet with those records only, in a new random order. The key is kept;
score uses the records present in the sheet. (Used on 2026-10-19 to cut the reading from 100 to 50 records.)

score: reads the filled sheet and the key and writes check/screening_check_results.txt and
deposit/screening_check.csv (check_id, DOI, title, the reader's decision and reason, the rules' decision).
Reported: the confusion matrix, Cohen's kappa on the sample, agreement in each stratum, and the overall agreement
weighted by the share of each stratum among all screened records, each with a 95% bootstrap interval
(resampling within strata, 10,000 replicates, seed 20261019).

Usage: run from data/ after pipeline/02_screen_records.py.
"""
import os
import sys

import numpy as np
import pandas as pd

DECISIONS_PATH = "pipeline_output/03_final_decisions.csv"
CHECK_DIR = "check"
SHEET = os.path.join(CHECK_DIR, "screening_check_sheet.xlsx")
KEY = os.path.join(CHECK_DIR, "screening_check_key.csv")
RESULTS = os.path.join(CHECK_DIR, "screening_check_results.txt")
DEPOSIT = "deposit/screening_check.csv"
SAMPLE_SIZE = {"Include": 50, "Maybe": 50, "Exclude": 35}
SHORT_SIZE = {"Include": 20, "Maybe": 15, "Exclude": 15}
SEED = 20261019
CLASSES = ["Include", "Maybe", "Exclude"]
REASONS = ["", "not_stochastic_framework", "not_genetics_domain", "keyword_coincidence", "wrong_document_type",
           "insufficient_information", "not_english", "duplicate_undetected", "other"]


def write_sheet(s):
    sheet = s[["check_id", "title", "abstract", "keywords", "source_title", "year", "doc_type"]].copy()
    sheet["reader_decision"] = ""
    sheet["reader_reason"] = ""
    sheet["reader_notes"] = ""
    with pd.ExcelWriter(SHEET, engine="openpyxl") as xw:
        pd.DataFrame({"Instructions": [
            "Read each record (title, abstract, keywords) and decide with the eligibility criteria below.",
            "reader_decision: Include, Maybe (title and abstract are not sufficient) or Exclude.",
            "reader_reason: for Exclude, one code from the list; leave empty for Include.",
            "Do not open screening_check_key.csv until every row is decided.",
            "",
            "Include: a stochastic model is developed or applied (I3), and genetics, evolution, plant or animal",
            "breeding, or gene-regulatory-network dynamics is the main biological subject (I4).",
            "Exclude: model applied outside these fields (not_genetics_domain); genetics mentioned only in passing",
            "(keyword_coincidence); deterministic method only (not_stochastic_framework); erratum, conference",
            "abstract or similar (wrong_document_type).",
        ]}).to_excel(xw, sheet_name="Instructions", index=False)
        sheet.to_excel(xw, sheet_name="Records", index=False)
        ws = xw.sheets["Records"]
        from openpyxl.worksheet.datavalidation import DataValidation
        for col, options in (("H", CLASSES), ("I", REASONS[1:])):
            dv = DataValidation(type="list", formula1='"' + ",".join(options) + '"', allow_blank=True)
            ws.add_data_validation(dv)
            dv.add(f"{col}2:{col}{len(sheet) + 1}")
        for col, width in zip("ABCDEFGHIJ", (8, 50, 90, 40, 30, 6, 10, 14, 24, 30)):
            ws.column_dimensions[col].width = width


def short():
    d = pd.read_csv(DECISIONS_PATH, encoding="utf-8-sig")
    key = pd.read_csv(KEY, encoding="utf-8-sig")
    pool = key[key.get("revealed_before_reading", pd.Series(0, index=key.index)).fillna(0).astype(int) == 0]
    rng = np.random.default_rng(SEED + 1)
    keep = []
    for cls in CLASSES:
        ids = pool.loc[pool["rules_decision"] == cls, "check_id"].to_numpy()
        keep.extend(rng.choice(ids, min(SHORT_SIZE[cls], len(ids)), replace=False))
    keep = rng.permutation(np.array(keep))
    rows = key.set_index("check_id").loc[keep, "record_row"]
    s = d.loc[rows.to_numpy()].reset_index(drop=True)
    s.insert(0, "check_id", keep)
    write_sheet(s)
    print(f"reading sheet shortened to {len(s)} records: "
          + ", ".join(f"{c} {int((s['final_decision'] == c).sum())}" for c in CLASSES))


def sample():
    d = pd.read_csv(DECISIONS_PATH, encoding="utf-8-sig")
    rng = np.random.default_rng(SEED)
    parts = []
    for cls in CLASSES:
        pool = d.index[d["final_decision"] == cls].to_numpy()
        n = min(SAMPLE_SIZE[cls], len(pool))
        parts.append(rng.choice(pool, n, replace=False))
    idx = rng.permutation(np.concatenate(parts))
    s = d.loc[idx].reset_index().rename(columns={"index": "record_row"})
    s.insert(0, "check_id", np.arange(1, len(s) + 1))
    os.makedirs(CHECK_DIR, exist_ok=True)
    s[["check_id", "record_row", "doi", "title", "year", "final_decision"]].rename(
        columns={"final_decision": "rules_decision"}).to_csv(KEY, index=False, encoding="utf-8-sig")

    write_sheet(s)
    print(f"sample of {len(s)} records: " + ", ".join(f"{c} {int((s['final_decision'] == c).sum())}" for c in CLASSES))
    print(f"wrote {SHEET} (to fill) and {KEY} (keep closed)")


def kappa(a, b):
    po = float(np.mean(a == b))
    pe = sum(float(np.mean(a == c)) * float(np.mean(b == c)) for c in CLASSES)
    return po, (po - pe) / (1 - pe) if pe < 1 else float("nan")


def score():
    key = pd.read_csv(KEY, encoding="utf-8-sig")
    filled = pd.read_excel(SHEET, sheet_name="Records")
    m = key.merge(filled[["check_id", "reader_decision", "reader_reason", "reader_notes"]], on="check_id")  # rows in the sheet
    # records whose rules' decision the reader saw before reading are not scored (column set by hand in the key)
    n_seen = int(m.get("revealed_before_reading", pd.Series(0, index=m.index)).fillna(0).astype(int).sum())
    if n_seen:
        m = m[m["revealed_before_reading"].fillna(0).astype(int) == 0].copy()
    m["reader_decision"] = m["reader_decision"].astype(str).str.strip().str.capitalize()
    missing = m[~m["reader_decision"].isin(CLASSES)]
    if len(missing):
        sys.exit(f"{len(missing)} rows have no valid reader_decision (check_id {missing['check_id'].tolist()[:10]} ...)")

    all_d = pd.read_csv(DECISIONS_PATH, encoding="utf-8-sig")["final_decision"]
    share = all_d.value_counts(normalize=True)
    a, b = m["reader_decision"].to_numpy(), m["rules_decision"].to_numpy()
    po, k = kappa(a, b)

    def weighted(df):
        return sum(share.get(c, 0) * float(np.mean(df.loc[df["rules_decision"] == c, "reader_decision"] == c))
                   for c in CLASSES if (df["rules_decision"] == c).any())

    rng = np.random.default_rng(SEED)
    strata = {c: m[m["rules_decision"] == c] for c in CLASSES}
    boot_k, boot_w = [], []
    for _ in range(10000):
        bs = pd.concat([g.iloc[rng.integers(0, len(g), len(g))] for g in strata.values() if len(g)])
        boot_k.append(kappa(bs["reader_decision"].to_numpy(), bs["rules_decision"].to_numpy())[1])
        boot_w.append(weighted(bs))
    ci = lambda x: (np.nanpercentile(x, 2.5), np.nanpercentile(x, 97.5))

    lines = [f"Screening check: {len(m)} records read blind to the rules' decision"
             + (f" ({n_seen} more left out: their rules' decision was seen before reading)" if n_seen else ""),
             "", "Confusion matrix (rows: reader, columns: rules)",
             pd.crosstab(m["reader_decision"], m["rules_decision"]).reindex(index=CLASSES, columns=CLASSES, fill_value=0).to_string(),
             "", f"Agreement in the sample: {po:.3f}; Cohen's kappa {k:.3f} (95% CI {ci(boot_k)[0]:.3f} to {ci(boot_k)[1]:.3f})"]
    for c in CLASSES:
        g = strata[c]
        if len(g):
            lines.append(f"Rules {c}: reader agreed on {int((g['reader_decision'] == c).sum())} of {len(g)}")
    w = weighted(m)
    lines += [f"Agreement weighted by stratum share among all {len(all_d)} records: {w:.3f} "
              f"(95% CI {ci(boot_w)[0]:.3f} to {ci(boot_w)[1]:.3f})"]
    open(RESULTS, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    m[["check_id", "doi", "title", "reader_decision", "reader_reason", "rules_decision"]].to_csv(
        DEPOSIT, index=False, encoding="utf-8-sig")
    print("\n".join(lines))
    print(f"\nwrote {RESULTS} and {DEPOSIT}")


if __name__ == "__main__":
    {"sample": sample, "short": short, "score": score}.get(sys.argv[1] if len(sys.argv) > 1 else "", lambda: sys.exit(__doc__))()
