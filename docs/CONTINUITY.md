# Course project status

The project is a standalone R collaboration for MCXA0046. Current scope is whole-grain exposure and longitudinal gut microbial variation in HMP2, with metabolomics available for context.

The data preparation entry point is `R/build_reference.R`. See `docs/DATA.md` for the input and pairing contract and `results/tables/reference_counts.csv` for the generated cohort summary.

Next analysis task: use the full metagenomic source to measure support for within-person changes in whole-grain responses and identify how zero profiles, repeated collections and time intervals affect eligibility. Then agree on the microbial outcome and sensitivity analysis before fitting a model.
