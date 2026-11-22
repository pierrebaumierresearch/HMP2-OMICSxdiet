# Whole-grain intake recency and gut microbiome change in Crohn’s disease

Does the gut microbiome of a person with Crohn’s disease change differently between stool samples when their report of whole-grain intake changes? This observational study uses repeated dietary questionnaires, stool metagenomes and stool metabolomics from the Inflammatory Bowel Disease Multi’omics Database (HMP2/IBDMDB). Each participant supplies their own comparison intervals, reducing the influence of stable differences between people.

## Study design

The primary analysis follows consecutive stool samples received 7–28 days apart. For each of 22 common genera, we compare two kinds of intervals from the same person that begin in the same dietary category. Switched intervals cross the “ate whole grains yesterday” boundary. Stable intervals stay on the same side. Reports of earlier intake and of no intake in the past seven days both fall outside “yesterday”. We estimate the two switch directions separately and give each participant equal weight. The exposure measures reported recency. The questionnaire does not measure portions or usual intake, and symptoms can influence food choice.

| Analysis | Question | Crohn’s cohort |
| --- | --- | --- |
| Whole grains and gut genera | How much does each genus change relative to the 22-genus panel on switched versus stable intervals? | 351 intervals, 35 participants |
| Wider dietary change | Is a larger number of changes across 21 dietary items associated with greater turnover of the 187-genus community? | 412 intervals, 52 participants |
| Stool butyrate | Within the same person, does stool butyrate intensity differ between visits when they report eating whole grains yesterday and visits when they report earlier or no intake in the past week? | 164 visits, 31 participants |
| Predictive benchmark | Does the updated whole-grain report improve estimation of the later microbial profile in a held-out participant, given the earlier profile, earlier report and time between samples? | 497 intervals, 52 participants |

The analysis plan specifies the genus-level contrasts, participant-level uncertainty intervals, Bayesian bootstrap and participant-held-out ridge benchmark. It also describes checks for changes in symptoms, medication and other dietary reports. [Statistical analysis plan](docs/HMP2-analysis-plan.md)

## Reproduce the data reconstruction

The R entry point builds a linked dataset for collections with both assays and complete responses for four foods and antibiotic use. It produces 341 collections from 76 participants across the Crohn’s and non-IBD groups. The analyses above use their own eligibility criteria and denominators. [Data sources and linkage](docs/DATA.md) · [Reconstruction counts](results/tables/reference_counts.csv)

From the repository root on Linux or WSL, use R and `sha256sum` to run the reconstruction.

```bash
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript R/build_reference.R --download
```

The build verifies the three input files against the SHA-256 hashes in the source manifest and writes the linked data to `data/processed/hmp2_reference.rds`. To use existing downloads, set `HMP2_RAW_DIR` to their directory and run `Rscript R/build_reference.R` without the download flag. With Pandoc installed, `Rscript -e 'rmarkdown::render("report/report.Rmd")'` renders the descriptive report. [Source manifest](data/sources.csv) · [R entry point](R/build_reference.R)

## Data citation

Lloyd-Price J, Arze C, Ananthakrishnan AN, et al. Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature* (2019). [Publication](https://doi.org/10.1038/s41586-019-1237-9) · [IBDMDB data downloads](https://ibdmdb.org/results)
