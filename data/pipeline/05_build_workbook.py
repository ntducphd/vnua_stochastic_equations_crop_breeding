"""
05_build_workbook.py -- PRISMA screening workbook (sheets Criteria, Screening, PRISMA_Summary,
Calibration_Validation).

The Decision / Exclusion_Reason / Notes / Reviewer columns hold the final decisions from
pipeline_output/03_final_decisions.csv. Every classifier suggestion was confirmed or overridden by the AI assistant
(Claude, Anthropic) at the first author's instruction, and the first author takes responsibility for the decisions
(manuscript Section 2.2). The classifier columns (orange header) are kept as an audit trail. Publisher abstracts
are not written to the workbook, because the database exports cannot be redistributed.

Usage: python pipeline/05_build_workbook.py   (run from data/, after step 04)
"""
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

IN_PATH = "pipeline_output/03_final_decisions.csv"
OUT_PATH = "pipeline_output/R01_PRISMA_screening_FINAL.xlsx"

df = pd.read_csv(IN_PATH, encoding="utf-8-sig")
n_rows = len(df)
n_changed = int((df["final_decision"] != df["ai_suggested_decision"]).sum())
n_exc_changed = int(((df["ai_suggested_decision"] == "Exclude") & (df["final_decision"] != "Exclude")).sum())
recon = pd.read_csv("pipeline_output/04_citation_reconciliation.csv", encoding="utf-8-sig")
n_cited, n_found = len(recon), int(recon["found_in_search"].sum())

wb = Workbook()

# ---------------- Sheet 1: Criteria ----------------
ws1 = wb.active
ws1.title = "Criteria"
ws1["A1"] = "Screening criteria (protocol/eligibility_criteria.md)"
ws1["A1"].font = Font(bold=True, size=13)
ws1.merge_cells("A1:B1")
crit_rows = [
    ("", ""),
    ("Inclusion criterion", ""),
    ("Development or application of a stochastic mathematical framework to genetics, "
     "evolution, plant or animal breeding, or gene regulatory network dynamics.", ""),
    ("", ""),
    ("Exclusion reasons (standardized labels used in Exclusion_Reason / ai_suggested_reason)", ""),
    ("not_stochastic_framework", "Deterministic methods only; no real stochastic/probabilistic framework"),
    ("not_genetics_domain", "Real stochastic framework, but applied to an unrelated domain "
                            "(linguistics/cultural evolution, computational genetic algorithms, econophysics, "
                            "ecology, epidemiology, etc.)"),
    ("keyword_coincidence", "Matched search terms only by incidental wording, a single passing mention, "
                            "or an ambiguous term with a non-genetics meaning (e.g. physical 'coalescence')"),
    ("duplicate_undetected", "A duplicate the DOI/title matching missed"),
    ("no_full_text", "Full text not accessible through any available channel"),
    ("not_english", "Full text not in English"),
    ("wrong_document_type", "Conference abstract, erratum, retraction, or similar non-substantive type"),
    ("insufficient_information", "Abstract missing or too sparse to judge"),
    ("other", "Any other reason - explain in Notes"),
]
r = 2
for a, b in crit_rows:
    ws1.cell(row=r, column=1, value=a)
    ws1.cell(row=r, column=2, value=b)
    if a and not b and r > 2:
        ws1.cell(row=r, column=1).font = Font(bold=True)
    r += 1
ws1.column_dimensions["A"].width = 45
ws1.column_dimensions["B"].width = 90
for row in ws1.iter_rows(min_row=2, max_row=r - 1):
    for cell in row:
        cell.alignment = Alignment(wrap_text=True, vertical="top")

# ---------------- Sheet 2: Screening (final) ----------------
ws2 = wb.create_sheet("Screening")
headers = ["ID", "Title", "Authors", "Year", "Source", "DOI", "Database", "Document_Type",
           "Abstract", "ai_suggested_decision", "ai_suggested_reason", "ai_confidence", "ai_flags",
           "Decision", "Exclusion_Reason", "Notes", "Reviewer"]
for col, h in enumerate(headers, start=1):
    c = ws2.cell(row=1, column=col, value=h)
    c.font = Font(bold=True, color="FFFFFF")
    is_ai_col = h.startswith("ai_")
    c.fill = PatternFill("solid", fgColor="ED7D31" if is_ai_col else "4472C4")
    c.alignment = Alignment(horizontal="center", vertical="center")

decision_fill = {
    "Include": PatternFill("solid", fgColor="E2EFDA"),
    "Exclude": PatternFill("solid", fgColor="FCE4E4"),
    "Maybe": PatternFill("solid", fgColor="FFF2CC"),
}

for i, row in df.iterrows():
    r = i + 2
    ws2.cell(row=r, column=1, value=i + 1)
    ws2.cell(row=r, column=2, value=row.get("title", ""))
    ws2.cell(row=r, column=3, value=row.get("authors", ""))
    ws2.cell(row=r, column=4, value=row.get("year", ""))
    ws2.cell(row=r, column=5, value=row.get("source_title", ""))
    ws2.cell(row=r, column=6, value=row.get("doi", ""))
    ws2.cell(row=r, column=7, value=row.get("database", ""))
    ws2.cell(row=r, column=8, value=row.get("doc_type", ""))
    ws2.cell(row=r, column=9, value="")   # publisher abstracts are not redistributed
    ws2.cell(row=r, column=10, value=row.get("ai_suggested_decision", ""))
    ws2.cell(row=r, column=11, value=row.get("ai_suggested_reason", ""))
    ws2.cell(row=r, column=12, value=row.get("ai_confidence", ""))
    ws2.cell(row=r, column=13, value=row.get("ai_flags", ""))
    for col in (10, 11, 12, 13):
        ws2.cell(row=r, column=col).font = Font(italic=True, color="7F7F7F")

    final_decision = row.get("final_decision", "")
    ws2.cell(row=r, column=14, value=final_decision)
    ws2.cell(row=r, column=15, value=row.get("final_exclusion_reason", ""))
    ws2.cell(row=r, column=16, value=row.get("final_notes", ""))
    ws2.cell(row=r, column=17, value="AI assistant (Claude, Anthropic); decisions accepted by Nguyen Trung Duc")
    if final_decision in decision_fill:
        ws2.cell(row=r, column=14).fill = decision_fill[final_decision]
    ws2.cell(row=r, column=14).font = Font(bold=True)

widths = [6, 42, 22, 6, 22, 20, 13, 11, 45, 16, 16, 10, 30, 12, 22, 45, 40]
for col, w in enumerate(widths, start=1):
    ws2.column_dimensions[get_column_letter(col)].width = w

ws2.freeze_panes = "B2"
ws2.auto_filter.ref = f"A1:Q{n_rows + 1}"

# ---------------- Sheet 3: PRISMA_Summary ----------------
ws3 = wb.create_sheet("PRISMA_Summary")
ws3["A1"] = "PRISMA Identification & Screening Summary - R01 (FINAL, rebuilt pipeline)"
ws3["A1"].font = Font(bold=True, size=13)
ws3.merge_cells("A1:C1")

ai_counts = df["ai_suggested_decision"].value_counts()
final_counts = df["final_decision"].value_counts()
exc_reasons = df[df["final_decision"] == "Exclude"]["final_exclusion_reason"].value_counts()

summary_rows = [
    ("Stage", "Count", "Source"),
    ("Records identified - Scopus", 974, "data/raw/scopus/ (export of 26 July 2026)"),
    ("Records identified - Web of Science", 721, "data/raw/WoS/ (export of 26 July 2026)"),
    ("Total records identified", 1695, "= Scopus + WoS"),
    ("Duplicate records removed", 1695 - n_rows, "deduplication by DOI, or by normalised title and year, "
                                        "pipeline/01_dedupe.py"),
    ("Unique records after deduplication", n_rows, "= Screening sheet row count"),
    ("", "", ""),
    ("AI-assisted screening (advisory only)", "", "pipeline/02_ai_screening.py, calibration kappa=0.365 (v2)"),
    ("  -> ai_suggested_decision = Include", int(ai_counts.get("Include", 0)), ""),
    ("  -> ai_suggested_decision = Maybe", int(ai_counts.get("Maybe", 0)), ""),
    ("  -> ai_suggested_decision = Exclude", int(ai_counts.get("Exclude", 0)), ""),
    ("", "", ""),
    ("FINAL screening decisions (classifier suggestions confirmed or overridden by the AI assistant)", "", "pipeline/03_merge_screening_decisions.py"),
    ("  -> Decision = Include", int(final_counts.get("Include", 0)), '=COUNTIF(Screening!N2:N' + str(n_rows + 1) + ',"Include")'),
    ("  -> Decision = Exclude", int(final_counts.get("Exclude", 0)), '=COUNTIF(Screening!N2:N' + str(n_rows + 1) + ',"Exclude")'),
    ("  -> Decision = Maybe (uncertain from title and abstract)", int(final_counts.get("Maybe", 0)),
     '=COUNTIF(Screening!N2:N' + str(n_rows + 1) + ',"Maybe")'),
    ("", "", ""),
    ("Exclude reasons breakdown", "", ""),
]
for reason, cnt in exc_reasons.items():
    summary_rows.append((f"  -> {reason}", int(cnt), ""))
summary_rows += [
    ("", "", ""),
    ("Records eligible for full-text assessment (Include + Maybe)",
     int(final_counts.get("Include", 0)) + int(final_counts.get("Maybe", 0)),
     "Next PRISMA stage - full-text eligibility review - NOT yet performed; "
     "this workbook covers title/abstract screening only"),
    ("Works cited in the manuscript", n_cited, "cited_references.csv"),
    ("  -> of which also retrieved by this Scopus/WoS search", n_found, "pipeline/04_reconcile_citations.py"),
    ("  -> of which not retrieved", n_cited - n_found,
     "mostly foundational theory or work outside the query terms; the search did not select the cited works"),
]
for i, (a, b, c) in enumerate(summary_rows):
    row = i + 3
    ws3.cell(row=row, column=1, value=a)
    ws3.cell(row=row, column=2, value=b)
    ws3.cell(row=row, column=3, value=c)
    if i == 0:
        for col in (1, 2, 3):
            ws3.cell(row=row, column=col).font = Font(bold=True)
            ws3.cell(row=row, column=col).fill = PatternFill("solid", fgColor="D9E1F2")
    if a.strip().startswith(("AI-assisted", "FINAL screening", "Exclude reasons")):
        ws3.cell(row=row, column=1).font = Font(bold=True, italic=True)

ws3.column_dimensions["A"].width = 60
ws3.column_dimensions["B"].width = 20
ws3.column_dimensions["C"].width = 65

# ---------------- Sheet 4: Calibration_Validation ----------------
ws4 = wb.create_sheet("Calibration_Validation")
ws4["A1"] = "AI Screening Calibration Validation"
ws4["A1"].font = Font(bold=True, size=13)
ws4.merge_cells("A1:D1")
cal_rows = [
    ("Metric", "v1 (initial)", "v2 (revised, used in this workbook)", "Note"),
    ("Validation sample size", 75, 75, "Same stratified sample (25 Include/25 Exclude/25 Maybe by v1 suggestion), "
                                       "read by the AI assistant without seeing the classifier suggestion"),
    ("Exact 3-way agreement", "45.3%", "70.7%", ""),
    ("Hard conflicts (Include vs Exclude flip)", 17, 2, "The metric that matters most for review integrity"),
    ("False Excludes of records the AI reader judged Include", 16, 0, "in-sample: v2 was tuned on these records"),
    ("Cohen's kappa", 0.180, 0.365, "'slight' -> 'fair' agreement; in-sample; the record-level readings were not kept"),
    ("Root cause diagnosed (v1)", "Required BOTH a domain-term AND a method-term match; several method terms "
     "(Wright-Fisher, chemical master equation, Moran model, coalescent, branching process, Potts model, "
     "diffusion approximation) are independently sufficient genetics-domain evidence in this corpus", "", ""),
    ("Residual known limitation (v2)", "Single incidental keyword mentions (e.g. a paper name-dropping "
     "'population genetics' once in a list of unrelated example domains) can still cause a false Include; "
     "not solved by keyword matching alone - addressed by the confirmation pass below", "", ""),
    ("", "", "", ""),
    ("CONFIRMATION PASS BY THE AI ASSISTANT (basis of the final decisions)", "", "", ""),
    ("Classifier suggestions confirmed or overridden by the AI assistant at the first author's instruction",
     n_rows, "every unique record; high-confidence inclusions checked from titles, in batches", ""),
    ("Decisions that differ from the classifier suggestion", n_changed, "documented with reason in pipeline/screening_decisions_ledger.csv, which also records the decisions on classifier-Maybe records"),
    ("Real example: high-confidence AI-Exclude bucket (n=15) still contained real errors",
     f"{n_exc_changed} of 15 changed", "e.g. a paper using 'genetic regulatory network' (not 'gene regulatory network') was "
     "wrongly excluded by the term list", ""),
]
for i, row_vals in enumerate(cal_rows):
    r = i + 3
    for col, v in enumerate(row_vals, start=1):
        c = ws4.cell(row=r, column=col, value=v)
        if i == 0:
            c.font = Font(bold=True, color="FFFFFF")
            c.fill = PatternFill("solid", fgColor="4472C4")
        if i == 9:
            c.font = Font(bold=True, italic=True)
        c.alignment = Alignment(wrap_text=True, vertical="top")
ws4.column_dimensions["A"].width = 42
ws4.column_dimensions["B"].width = 45
ws4.column_dimensions["C"].width = 45
ws4.column_dimensions["D"].width = 45

wb.save(OUT_PATH)
print(f"Saved {OUT_PATH} - {n_rows} records, final decisions (confirmed by the AI assistant) populated.")
print(final_counts)
