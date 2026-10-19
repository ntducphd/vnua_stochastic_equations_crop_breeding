# Screening check: analysis plan

Fixed on 19 October 2026, before the key of the check sample was opened. The commit that adds this file to the
code repository dates it.

## Design

- **Rules under test:** `data/pipeline/screening_rules.yaml`, version 1.2. The rules were frozen after being
  refined on a random half of the records.
- **Sample:** a random sample stratified by the rules' decision (seed 20261019), shortened with seed 20261020 to
  20 Include, 15 Maybe and 15 Exclude records. The scripts are `data/pipeline/03_check_screening.py sample` and
  `short`. Records whose rules' decision the reader saw before reading are left out of the scoring.
- **Reader:** the first author, who reads title, abstract and keywords only, blind to the rules' decisions, and
  applies `eligibility_criteria.md`. The decision is Include, Maybe or Exclude, with an exclusion reason.
- There is one reader and no second, independent reader. This is a limitation of the check.

## Metrics

Estimates refer to all 1,059 screened records. Each stratum is weighted by its share of the screened records.

- **Primary:** recall, the probability that the rules do not exclude a record the reader includes.
- **Acceptance:** the rules are accepted when the recall estimate is at least **0.95**.
- **Secondary:**
  - strict recall (the rules include the record);
  - specificity (the rules exclude a record the reader excludes);
  - agreement over the three classes;
  - the estimated number of eligible records among those the rules exclude;
  - Cohen's kappa and PABAK on the sample.
- Every estimate carries a 95% bootstrap percentile interval: 10,000 replicates, resampled within strata,
  seed 20261019.
- **Reader Maybe:** counts in agreement and kappa, but not in the denominator of recall or of specificity.

## Reporting

- The paper reports every metric above, including any that falls short of the threshold.
- If recall is below 0.95, the paper says so. Any change to the rules made after the key is opened is then
  reported as a post hoc revision, and the check is not re-scored on the same sample as if it were independent.
