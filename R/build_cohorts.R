source("R/01_sources.R")
source("R/02_cohort.R")
source("R/03_metabolomics.R")
source("R/04_longitudinal.R")
source("R/05_butyrate.R")

paths <- study_inputs("R/build_cohorts.R")
metadata_path <- paths[["hmp2_metadata_2018-08-20.csv"]]
mgx <- build_mgx_cohorts(metadata_path, paths[["taxonomic_profiles_3.tsv.gz"]])
mbx <- build_butyrate_cohort(metadata_path, paths[["HMP2_metabolomics_w_metadata.biom.gz"]])
flow <- rbind(mgx$flow, mbx$flow)

# These four denominators were specified using the pinned source releases.
expected <- data.frame(
  stage = c("mgx_primary", "mgx_whole_recall", "mbx_butyrate_both_states", "mgx_all_intervals"),
  unit = c("interval", "interval", "visit", "interval"),
  n = c(351L, 412L, 164L, 497L),
  participants = c(35L, 52L, 31L, 52L)
)
observed <- flow[match(expected$stage, flow$stage), names(expected), drop = FALSE]
if (anyDuplicated(flow$stage) || anyNA(observed) || !all(observed == expected)) {
  stop("Analytic cohort denominators differ from the protocol. No files written.")
}

members <- list(mgx_intervals = mgx$intervals, genus_panel = mgx$panel,
                butyrate_visits = mbx$visits)
write_output("data/processed/study_cohorts.rds", function(path) saveRDS(members, path))
write_output("results/tables/cohort_flow.csv", function(path) {
  write.csv(flow, path, row.names = FALSE)
})
print(observed, row.names = FALSE)
