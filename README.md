# Stochastic equations in crop breeding

Code and data for the review

> Nguyen Trung Duc, Nguyen Thi Hoai Linh, Ta Viet Ton. Stochastic equations in crop breeding: a review of
> models and applications. Manuscript in preparation.

## Contents

| Folder | Contents | In the paper |
|---|---|---|
| `scripts/` | Wright–Fisher, diffusion and branching computations; model figures | Figures 1, 4 and 5; Table 2; Box 2; Sections 5 and 7 |
| `data/` | literature search: pipeline, screening decisions, evidence-map tables | Figures 2 and 3; Section 2 |
| `genome/` | *a1*–*sh2* distances and chromosome 4 diversity | Sections 3.3 and 4.4 |
| `protocol/` | eligibility criteria of the literature search | Section 2 |
| `figures/` | figures as PDF and 600 dpi PNG | |

## Requirements

Python 3.13 and the packages in `requirements.txt`:

    pip install -r requirements.txt

## Running

From the repository root:

| Command | Output |
|---|---|
| `python -I -B scripts/generate_figures.py` | Figures 1, 4 and 5, and the values quoted in their captions |
| `python -I -B scripts/wf_exact.py` | exact Wright–Fisher values of Table 2 and Box 2, and the *a1*–*sh2* screening sizes |
| `python -I -B scripts/worked_examples.py` | diffusion times, branching probabilities and replicate numbers; writes `scripts/worked_examples_output.txt` |
| `python -I -B scripts/monte_carlo_check.py` | Monte Carlo check of Table 2 (100,000 replicates per scenario by default) |
| `python -I -B scripts/check_eigenvalue_inflation.py` | largest-eigenvalue inflation of an estimated covariance matrix (Table 3) |
| `python -I -B genome/scripts/genome_numbers.py` | the genome values quoted in Sections 3.3 and 4.4 |
| `cd data && python -I -B pipeline/test_screen_records.py` | one test case for each screening rule (needs no database export) |

The Euler–Maruyama paths of Figure 5a use seed 42; the other figures are deterministic. The literature-search
pipeline is described in `data/README.md` and the genome tables in `genome/README.md`.

## Not included

The Scopus and Web of Science exports and the files derived from them contain publisher abstracts and author
contact details and cannot be redistributed. `data/deposit/records_identified.csv` lists the identifiers from
which they can be rebuilt.

## Licence

Code: MIT (`LICENSE`). Tables created for this study (`data/deposit/`, `data/pipeline/*.csv`,
`genome/data/curated/`, `genome/results/`): CC BY 4.0. Bibliographic data of the screened records come from
Scopus and Web of Science; the genome tables are derived from the data of Hu et al. (2021, CC BY) and MaizeGDB.
