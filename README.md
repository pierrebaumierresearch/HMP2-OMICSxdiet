# Diet and gut profiles in Crohn's disease

Collaborative mini-project for MCXA0046, AI in Health Environment. We use the public HMP2/IBDMDB study to examine whether reported recent whole-grain consumption is associated with within-person variation in gut microbial profiles. The metabolome is available for biological context. The other dietary questions and the antibiotic flag remain available for secondary checks.

## Reproduce the reference table

The current pipeline provides descriptive data preparation. The statistical analysis and written report will follow the course protocol once the team defines the analysis population and outcome.

From the project root, on Linux or WSL with R and `sha256sum` available:

```bash
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript R/build_reference.R --download
```

With Pandoc available, run `Rscript -e 'rmarkdown::render("report/report.Rmd")'` to generate the Markdown report.

R reads the three public products listed in `data/sources.csv`. Use `HMP2_RAW_DIR` to point to an existing folder containing all three files under their original names. The script verifies their SHA-256 hashes. A run without `--download` uses existing inputs only.

| Location | Contents |
| --- | --- |
| `R/` | Source verification, exposure and cohort construction, BIOM parser and build entry point |
| `data/raw/` | Local source downloads |
| `data/processed/hmp2_reference.rds` | Collection metadata, genus matrix, peak matrix and peak annotations |
| `results/tables/` | Generated summary tables |
| `results/figures/` | Analysis figures |
| `report/` | Written course report in R Markdown |
| `docs/` | Provenance and project status |

The RDS stores each omic matrix by collection order alongside the selected assay identifiers. Git tracks the aggregate cohort counts in `results/tables/reference_counts.csv`. Input files and generated participant-level datasets stay in local ignored directories. The [IBDMDB download page](https://ibdmdb.org/results) provides the source data and metadata.
