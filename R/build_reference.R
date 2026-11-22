source("R/01_sources.R")
source("R/02_cohort.R")
source("R/03_metabolomics.R")

paths <- study_inputs("R/build_reference.R")

cohort <- build_reference_cohort(paths[["hmp2_metadata_2018-08-20.csv"]])
mgx <- read_genus(paths[["taxonomic_profiles_3.tsv.gz"]], cohort$mgx_sample_id)
mbx <- read_metabolomics(paths[["HMP2_metabolomics_w_metadata.biom.gz"]], cohort$mbx_sample_id)
cohort$microbiome_all_zero <- rowSums(mgx != 0) == 0
cohort$metabolome_all_zero <- rowSums(mbx$intensity != 0) == 0
stopifnot(identical(rownames(mgx), cohort$mgx_sample_id),
          identical(rownames(mbx$intensity), cohort$mbx_sample_id))

people <- !duplicated(cohort$participant_id)
observations <- table(cohort$participant_id)
summary <- c(
  collections = nrow(cohort), participants = sum(people),
  crohn_participants = sum(cohort$diagnosis[people] == "CD"),
  non_ibd_participants = sum(cohort$diagnosis[people] == "nonIBD"),
  crohn_collections = sum(cohort$diagnosis == "CD"),
  non_ibd_collections = sum(cohort$diagnosis == "nonIBD"),
  participants_with_repeat_collections = sum(observations >= 2L),
  exact_id_pairs = sum(cohort$pairing_rule == "exact_external_id"),
  documentary_collection_pairs = sum(cohort$pairing_rule == "unique_documentary_collection_key"),
  microbiome_genera = ncol(mgx), metabolomic_peaks = ncol(mbx$intensity),
  microbiome_all_zero = sum(cohort$microbiome_all_zero),
  metabolome_all_zero = sum(cohort$metabolome_all_zero)
)
expected <- c(341L, 76L, 50L, 26L, 218L, 123L, 73L, 279L, 62L, 187L, 81867L, 2L, 0L)
if (!identical(as.integer(summary), expected)) stop("Reference cohort or matrix dimensions changed")

write_output("data/processed/hmp2_reference.rds", function(path) {
  saveRDS(list(cohort = cohort, genus = mgx, metabolome = mbx$intensity,
               metabolome_features = mbx$features), path, compress = "gzip")
})
write_output("results/tables/reference_counts.csv", function(path) {
  write.csv(data.frame(metric = names(summary), value = unname(summary)), path, row.names = FALSE)
})
print(summary)
