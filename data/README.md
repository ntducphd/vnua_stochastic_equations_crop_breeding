# Literature search

Scopus and Web of Science Core Collection were searched without a lower date limit and exported on 26 July
2026: 1,695 records (974 Scopus, 721 Web of Science), 1,059 after deduplication. The search terms are in
Section 2.1 of the paper; the exact query strings were not kept. Eligibility criteria:
`../protocol/eligibility_criteria.md`.

Titles, abstracts and keywords are screened by explicit rules (step 02). `pipeline/screening_rules.yaml`
writes each eligibility criterion as term lists and thresholds; `pipeline/02_screen_records.py` applies them and
records, for every record, the decision, the rule that made it, the matched terms and the SHA-256 of the rules
file. The same records and rules always give the same decisions. The rules were refined on a random half of the
records and then frozen (version 1.2). Step 03 checks them against a blinded human reading of a stratified random
sample of 100 records. Full texts were not assessed.

## Pipeline

Run from this folder, in order:

| Step | Script | Output |
|---|---|---|
| 01 | `pipeline/01_dedupe.py` | deduplicated records and identification counts |
| 02 | `pipeline/02_screen_records.py` | screening decisions, rules and matched terms (`pipeline/screening_rules.yaml`) |
| 03 | `pipeline/03_check_screening.py sample` / `score` | blinded check sample of 100 records, then its agreement with the rules |
| test | `pipeline/test_screen_records.py` | one test case per screening rule |
| 04 | `pipeline/04_reconcile_citations.py` | cited works also retrieved by the search |
| 05 | `pipeline/05_build_workbook.py` | screening workbook |
| 06 | `pipeline/06_render_figure1_prisma_flow.py` | Figure 2 (PRISMA 2020 flow diagram) |
| 12 | `pipeline/12_build_deposit_tables.py` | `deposit/` |
| 13 | `pipeline/13_evidence_map.py` | Figure 3 (evidence map) and its counts |

Every step reads the database exports (`raw/`) or files derived from them. The exports contain publisher
abstracts and author contact details and are not redistributed; they can be rebuilt from the identifiers in
`deposit/records_identified.csv`.

Records are deduplicated by matching DOI, or matching title (letters and digits only) and year. Step 13
assigns equation families and application domains with the regular expressions at the top of the script.

## Deposit

| File | Contents |
|---|---|
| `deposit/records_identified.csv` | the 1,695 retrieved records: database, Scopus EID or Web of Science UT, DOI, title, year, source, document type |
| `deposit/screening_decisions.csv` | the 1,059 deduplicated records with the decision, reason, rule, scores, matched terms and rules version |
| `deposit/screening_check.csv` | the blinded reading of the check sample next to the rules' decision (written by step 03 `score`) |
| `deposit/citation_reconciliation.csv` | the 135 cited works matched against the screened records |
| `deposit/R01_screening_workbook.xlsx` | screening workbook: criteria, decisions, PRISMA counts, calibration summary |
| `deposit/evidence_map_records.csv` | family and domain flags of the 989 retained records, with identifiers and year |
| `deposit/evidence_map_counts.csv` | counts by family, domain and period (long format) |
| `cited_references.csv` | the reference list of the paper, read by step 04 |

## Results

- Screening decisions (rules 1.2): 1,000 retained, 15 uncertain, 44 excluded (37 outside genetics or breeding,
  4 matched through another sense of the words, 3 without a stochastic model).
- 9 of the 135 cited works were also retrieved by the search.
- Records per family: master equations and Markov chains 279, diffusion 680, coalescent 168, branching 202,
  control 8, simulation 79, no family term 43. Plant terms occur in 46 records and breeding terms in 15.
