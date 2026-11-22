# Whole-grain intake recency and gut microbiome change in Crohn’s disease

Does the gut microbiome of a person with Crohn’s disease change differently between stool samples when their report of whole-grain intake changes? This observational study uses repeated dietary questionnaires, stool metagenomes and stool metabolomics from the Inflammatory Bowel Disease Multi’omics Database (HMP2/IBDMDB). Each participant supplies their own comparison intervals, reducing the influence of stable differences between people.

## Study design

The primary analysis follows consecutive stool samples received 7–28 days apart. For each of 22 common genera, we compare two kinds of intervals from the same person that begin in the same dietary category. Switched intervals cross the “ate whole grains yesterday” boundary. Stable intervals stay on the same side. Reports of earlier intake and of no intake in the past seven days both fall outside “yesterday”. We estimate the two switch directions separately and give each participant equal weight. The exposure measures reported recency. The questionnaire does not measure portions or usual intake, and symptoms can influence food choice.

The R scripts reproduce eligibility and membership for all four cohorts below. Association estimates and predictive performance are specified in the protocol. [Cohort flow](results/tables/cohort_flow.csv)

| Analysis | Question | Crohn’s cohort |
| --- | --- | --- |
| Whole grains and gut genera | How much does each genus change relative to the 22-genus panel on switched versus stable intervals? | 351 intervals, 35 participants |
| Wider dietary change | Is a larger number of changes across 21 dietary items associated with greater turnover of the 187-genus community? | 412 intervals, 52 participants |
| Stool butyrate | Within the same person, does stool butyrate intensity differ between visits when they report eating whole grains yesterday and visits when they report earlier or no intake in the past week? | 164 visits, 31 participants |
| Predictive benchmark | Does the updated whole-grain report improve estimation of the later microbial profile in a held-out participant, given the earlier profile, earlier report and time between samples? | 497 intervals, 52 participants |

The analysis plan specifies the genus-level contrasts, participant-level uncertainty intervals, Bayesian bootstrap and participant-held-out ridge benchmark. It also describes checks for changes in symptoms, medication and other dietary reports. [Statistical analysis plan](docs/HMP2-analysis-plan.md)

## Reproduce the study cohorts

From the repository root on Linux or WSL, use R and `sha256sum`. The lockfile records R 4.5.0 and the package versions.

```bash
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript R/build_cohorts.R --download
```

This command verifies the three input files, writes the aggregate counts to `results/tables/cohort_flow.csv`, and stores participant-level interval and visit membership locally in `data/processed/study_cohorts.rds`. To reuse existing downloads, set `HMP2_RAW_DIR` to their directory and omit `--download`. [Source manifest](data/sources.csv) · [Dietary item definitions](data/diet_items.csv) · [R entry point](R/build_cohorts.R)

A separate reconstruction describes 341 collections from 76 participants across the Crohn’s and non-IBD groups with both assays, consistent answers for four foods and a consistent exported antibiotic flag. A blank parent medication question makes a negative flag ambiguous. With the inputs available, run `Rscript R/build_reference.R` to reproduce this table. With Pandoc installed, `Rscript -e 'rmarkdown::render("report/report.Rmd")'` renders its descriptive report. [Linkage criteria](docs/DATA.md) · [Linkage counts](results/tables/reference_counts.csv)

## Contributing

Use a short branch and open a pull request for a change to the analysis. Include the command you ran and any change in aggregate counts. Run `Rscript --vanilla tests/test_study_cohorts.R` to check the cohort rules on synthetic data. Keep participant-level files local and share aggregate tables through `results/tables/`.

## Data citation

Lloyd-Price J, Arze C, Ananthakrishnan AN, et al. Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature* (2019). [Publication](https://doi.org/10.1038/s41586-019-1237-9) · [IBDMDB data downloads](https://ibdmdb.org/results)
