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
The metrics and the acceptance threshold were fixed before the key was opened (../protocol/screening_check_plan.md).
Population quantities are estimated with stratum weights: each stratum (the rules' decision) counts with its share
of all screened records, so the over-sampled Maybe and Exclude strata do not bias the estimates.
  recall (primary)    P(rules do not exclude | reader Include); accepted when the estimate is >= RECALL_THRESHOLD
  strict recall       P(rules Include | reader Include)
  specificity         P(rules Exclude | reader Exclude)
  weighted agreement  P(rules decision = reader decision), three classes
  missed eligible     estimated number of reader-Include records among all records the rules exclude
  Cohen's kappa and PABAK (prevalence- and bias-adjusted kappa, three classes) on the sample
Reader Maybe counts in the agreement and kappa but in neither the recall nor the specificity denominator. Every
interval is a 95% bootstrap percentile interval (resampling within strata, 10,000 replicates, seed 20261019).

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
RECALL_THRESHOLD = 0.95          # pre-specified, ../protocol/screening_check_plan.md
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


def pabak(a, b):
    k = len(CLASSES)
    return (k * float(np.mean(a == b)) - 1) / (k - 1)


def weighted_metrics(df, share, n_rules):
    """Population estimates from a sample stratified by the rules' decision."""
    joint = {}                                    # P(rules = s, reader = r)
    for s_ in CLASSES:
        g = df[df["rules_decision"] == s_]
        for r in CLASSES:
            joint[s_, r] = share.get(s_, 0) * (float(np.mean(g["reader_decision"] == r)) if len(g) else 0.0)
    elig = sum(joint[s_, "Include"] for s_ in CLASSES)
    inel = sum(joint[s_, "Exclude"] for s_ in CLASSES)
    nan = float("nan")
    excl = df[df["rules_decision"] == "Exclude"]
    return {
        "recall": (joint["Include", "Include"] + joint["Maybe", "Include"]) / elig if elig else nan,
        "strict_recall": joint["Include", "Include"] / elig if elig else nan,
        "specificity": joint["Exclude", "Exclude"] / inel if inel else nan,
        "agreement": sum(joint[c, c] for c in CLASSES),
        "missed": n_rules.get("Exclude", 0) * (float(np.mean(excl["reader_decision"] == "Include")) if len(excl) else 0.0),
    }


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
    share, n_rules = all_d.value_counts(normalize=True), all_d.value_counts()
    a, b = m["reader_decision"].to_numpy(), m["rules_decision"].to_numpy()
    po, k = kappa(a, b)
    pb = pabak(a, b)
    est = weighted_metrics(m, share, n_rules)

    rng = np.random.default_rng(SEED)
    strata = [m[m["rules_decision"] == c] for c in CLASSES if (m["rules_decision"] == c).any()]
    boot = {name: [] for name in list(est) + ["kappa", "pabak"]}
    for _ in range(10000):
        bs = pd.concat([g.iloc[rng.integers(0, len(g), len(g))] for g in strata])
        x, y = bs["reader_decision"].to_numpy(), bs["rules_decision"].to_numpy()
        boot["kappa"].append(kappa(x, y)[1])
        boot["pabak"].append(pabak(x, y))
        for name, v in weighted_metrics(bs, share, n_rules).items():
            boot[name].append(v)
    ci = lambda name: (np.nanpercentile(boot[name], 2.5), np.nanpercentile(boot[name], 97.5))
    fmt = lambda name, v, d=3: f"{v:.{d}f} (95% CI {ci(name)[0]:.{d}f} to {ci(name)[1]:.{d}f})"

    lines = [f"Screening check: {len(m)} records read blind to the rules' decision"
             + (f" ({n_seen} more left out: their rules' decision was seen before reading)" if n_seen else ""),
             "", "Confusion matrix (rows: reader, columns: rules)",
             pd.crosstab(m["reader_decision"], m["rules_decision"]).reindex(index=CLASSES, columns=CLASSES, fill_value=0).to_string(),
             ""]
    for c in CLASSES:
        g = m[m["rules_decision"] == c]
        if len(g):
            lines.append(f"Rules {c}: reader agreed on {int((g['reader_decision'] == c).sum())} of {len(g)}")
    lines += ["", f"Estimates for all {len(all_d)} screened records (stratum-weighted):",
              f"  recall (rules do not exclude a reader-Include record): {fmt('recall', est['recall'])}"
              f"  -> {'meets' if est['recall'] >= RECALL_THRESHOLD else 'does NOT meet'} the pre-specified threshold {RECALL_THRESHOLD}",
              f"  strict recall (rules Include):                        {fmt('strict_recall', est['strict_recall'])}",
              f"  specificity (rules Exclude a reader-Exclude record):   {fmt('specificity', est['specificity'])}",
              f"  agreement, three classes:                              {fmt('agreement', est['agreement'])}",
              f"  eligible records among the {int(n_rules.get('Exclude', 0))} the rules exclude: "
              f"{fmt('missed', est['missed'], 1)}",
              "", f"On the sample: agreement {po:.3f}; Cohen's kappa {fmt('kappa', k)}; PABAK {fmt('pabak', pb)}"]
    open(RESULTS, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    m[["check_id", "doi", "title", "reader_decision", "reader_reason", "rules_decision"]].to_csv(
        DEPOSIT, index=False, encoding="utf-8-sig")
    print("\n".join(lines))
    print(f"\nwrote {RESULTS} and {DEPOSIT}")


if __name__ == "__main__":
    {"sample": sample, "short": short, "score": score}.get(sys.argv[1] if len(sys.argv) > 1 else "", lambda: sys.exit(__doc__))()
