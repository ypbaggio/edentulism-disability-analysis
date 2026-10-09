# Analysis Code for Edentulism Disability and Death Across Six International Aging Cohorts

This repository contains the statistical analysis code and deidentified aggregate results supporting the manuscript *Edentulism, Disability, and Death Across Six International Aging Cohorts*.

## Scope

The public workflow starts from cohort-specific harmonized first-event analysis panels. It includes the primary multiple-imputation models, secondary one-stage models, sensitivity models, and the six-cohort random-effects synthesis. Raw cohort files, participant-level analysis files, and raw-data ingestion or harmonization programs are not distributed.

CHARLS, CLHLS, KLoSA, HRS, ELSA, and SHARE remain subject to their respective registration, access, and data-use terms. Researchers seeking to rerun participant-level models must obtain the source data independently and construct a panel that follows `data/input_schema.csv` and the definitions in the manuscript and Supplemental Appendix.

No participant identifiers or participant-level records are included in this repository.

## Reproducible aggregate analysis

The primary six-cohort synthesis can be reproduced without restricted data:

```bash
python analysis/04_pool_primary_estimates.py --check
```

The script reads `data/primary_cohort_estimates.csv`, performs restricted-maximum-likelihood random-effects meta-analysis with modified Hartung-Knapp inference, calculates 95% prediction intervals, writes `output/reproduced_primary_pooled_results.csv`, and verifies the result against the reported estimates.

## Participant-level analysis code

The Stata programs require a locally held harmonized first-event panel and never write participant-level data to this repository.

1. `analysis/01_primary_mi_models.do` performs 20-fold chained multiple imputation separately by cohort and fits the cohort-specific Poisson models.
2. `analysis/02_one_stage_models.do` fits the secondary mixed-effects models with cohort-specific attained-age functions and t inference with 5 degrees of freedom.
3. `analysis/03_one_stage_sensitivities.do` fits fixed-effects, observation-weighted, complete-case, and random-slope sensitivity models.

Example Stata usage:

```stata
do "analysis/01_primary_mi_models.do" "LOCAL_ANALYSIS_PANEL.dta" "LOCAL_OUTPUT_DIRECTORY"
do "analysis/02_one_stage_models.do" "LOCAL_ANALYSIS_PANEL.dta" "LOCAL_OUTPUT_DIRECTORY"
do "analysis/03_one_stage_sensitivities.do" "LOCAL_ANALYSIS_PANEL.dta" "LOCAL_OUTPUT_DIRECTORY"
```

## Repository contents

- `analysis/`: statistical model and synthesis code.
- `data/input_schema.csv`: required structure of the harmonized first-event panel.
- `data/primary_cohort_estimates.csv`: cohort-specific multiply imputed estimates used in the primary synthesis.
- `data/reported_primary_pooled_results.csv`: reported six-cohort pooled estimates.
- `data/reported_*` and selected descriptive CSV files: deidentified aggregate results used to check major sensitivity and secondary analyses.
- `output/`: generated results; the directory is empty before reproduction.

## Software

- Python 3.12 or later with NumPy, pandas, and SciPy.
- Stata 19 for participant-level models.

Install the Python dependencies with:

```bash
python -m pip install -r requirements.txt
```

## Interpretation boundary

The repository enables direct verification of the cross-cohort synthesis and inspection of the participant-level model specifications. Because the source cohort data cannot be redistributed, it does not by itself recreate the harmonized participant-level panels from raw custodian files.

