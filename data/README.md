# Literature search

Scopus and Web of Science Core Collection were searched without a lower date limit and exported on 26 July
2026: 1,695 records (974 Scopus, 721 Web of Science), 1,059 after deduplication. The search terms are in
Section 2.1 of the paper; the exact query strings were not kept. Eligibility criteria:
`../protocol/eligibility_criteria.md`.

A rule-based classifier (step 02) suggested a decision for each record. An AI assistant (Claude, Anthropic)
confirmed or overrode each suggestion at the first author's instruction and changed 177 of them. On 75
records read by the assistant without the suggestions, the classifier agreed with it at κ = 0.365 after
revision (0.180 before); both values are in-sample. Full texts were not assessed.

## Pipeline

Run from this folder, in order:

| Step | Script | Output |
|---|---|---|
| 01 | `pipeline/01_dedupe.py` | deduplicated records and identification counts |
| 02 | `pipeline/02_ai_screening.py` | classifier suggestions |
| 03 | `pipeline/03_merge_screening_decisions.py` | final decisions (applies `pipeline/screening_decisions_ledger.csv`) |
| 04 | `pipeline/04_reconcile_citations.py` | cited works also retrieved by the search |
| 05 | `pipeline/05_build_workbook.py` | screening workbook |
| 06 | `pipeline/06_render_figure1_prisma_flow.py` | Figure 2 (PRISMA 2020 flow diagram) |
| 12 | `pipeline/12_build_deposit_tables.py` | `deposit/` |
| 13 | `pipeline/13_evidence_map.py` | Figure 3 (evidence map) and its counts |

Every step reads the database exports (`raw/`) or files derived from them. The exports contain publisher
abstracts and author contact details and are not redistributed; they can be rebuilt from the identifiers in
`deposit/records_identified.csv`.

Records are deduplicated by matching DOI, or matching title (letters and digits only) and year. Step 13
assigns equation families and application domains with the regular expressions at the top of the script;
`pipeline/13_evidence_map_handcheck.csv` holds the hand labels of 60 records used to check them.

## Deposit

| File | Contents |
|---|---|
| `deposit/records_identified.csv` | the 1,695 retrieved records: database, Scopus EID or Web of Science UT, DOI, title, year, source, document type |
| `deposit/screening_decisions.csv` | the 1,059 deduplicated records with the classifier suggestion and the final decision and reason |
| `deposit/screening_decisions_ledger.csv` | the 202 decisions entered through the ledger, with reasons |
| `deposit/citation_reconciliation.csv` | the 135 cited works matched against the screened records |
| `deposit/R01_screening_workbook.xlsx` | screening workbook: criteria, decisions, PRISMA counts, calibration summary |
| `deposit/evidence_map_records.csv` | family and domain flags of the 989 retained records, with identifiers and year |
| `deposit/evidence_map_counts.csv` | counts by family, domain and period (long format) |
| `cited_references.csv` | the reference list of the paper, read by step 04 |

## Results

- Final decisions: 989 retained, 27 uncertain, 43 excluded (30 outside genetics or breeding, 8 incidental
  keyword matches, 5 without a stochastic model).
- 9 of the 135 cited works were also retrieved by the search.
- Records per family: master equations and Markov chains 282, diffusion 675, coalescent 167, branching 196,
  control 9, simulation 78, no family term 40. Plant terms occur in 42 records and breeding terms in 13.
