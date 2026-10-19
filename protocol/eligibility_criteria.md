# Eligibility criteria

Title-and-abstract screening of the records retrieved from Scopus and Web of Science Core Collection
(exported 26 July 2026) for the review *Stochastic equations in crop breeding: a review of models and
applications*.

## Inclusion

1. Peer-reviewed journal article, review article or preprint with a DOI.
2. Written in English.
3. Develops or applies at least one stochastic model, for example a chemical master equation, a
   Wright–Fisher chain or diffusion, a Fokker–Planck equation, a branching process, Hamilton–Jacobi–Bellman
   optimal control, a mean-field game, an Ising or Potts analogy for genetics, or a neural or score-based SDE.
4. Genetics, evolution, plant or animal breeding, or gene-regulatory-network dynamics is the main
   biological subject.
5. Title and abstract are sufficient for a decision.

## Exclusion

1. The stochastic model is applied outside genetics, evolution, breeding and gene regulation (finance,
   epidemiology without population genetics, image processing, robotics, quantum physics, engineering control).
2. Genetics or evolution is mentioned only in passing.
3. The method is deterministic.
4. Duplicate record.
5. Not in English, or no usable abstract.
6. Conference abstract, erratum, retraction notice or another non-substantive document.

## Exclusion reasons

| Code | Meaning |
|---|---|
| `not_stochastic_framework` | Deterministic methods only |
| `not_genetics_domain` | Stochastic model outside the biological scope |
| `keyword_coincidence` | Search terms matched by incidental wording |
| `duplicate_undetected` | Duplicate missed by DOI and title matching |
| `no_full_text` | Text not accessible |
| `not_english` | Not in English |
| `wrong_document_type` | Conference abstract, erratum or similar |
| `insufficient_information` | Abstract missing or too short |
| `other` | Explained in the notes column |

## Screening

The criteria are applied by explicit rules: `data/pipeline/screening_rules.yaml` holds the term lists and
thresholds for criteria I3, I4, I5, E1, E2, E3, E5 and E6, and `data/pipeline/02_screen_records.py` gives each
record a decision, a reason, the rule that decided it and the matched terms. Records with too little evidence
either way are classed as uncertain (Maybe). A stratified random sample of 100 records is read blind to the rules'
decisions (`data/pipeline/03_check_screening.py`). Full texts were not assessed. Decisions are in
`data/deposit/screening_decisions.csv` and the check in `data/deposit/screening_check.csv`. The analysis plan of the check, with its
acceptance threshold, was fixed before the key was opened: `screening_check_plan.md`.
