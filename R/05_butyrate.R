butyrate_membership <- function(visits, intensity) {
  # Join the peak by assay ID. A zero in this sparse export is unusable on the
  # log scale; an absent assay column is a linkage error.
  keys <- visits[c("participant_id", "sample_id", "collection_id")]
  stopifnot(!anyNA(keys), all(vapply(keys, function(x) all(nzchar(x)), logical(1))),
            !anyDuplicated(visits$sample_id),
            !anyDuplicated(visits[c("participant_id", "collection_id")]),
            is.numeric(intensity), !is.null(names(intensity)), !anyDuplicated(names(intensity)))
  matched <- match(visits$sample_id, names(intensity))
  stopifnot(!anyNA(matched))
  values <- intensity[matched]
  positive <- is.finite(values) & values > 0
  answered <- positive & visits$whole_grain %in% food_answers
  visits$yesterday <- visits$whole_grain %in% tail(food_answers, 2L)
  people_y <- visits$participant_id[answered & visits$yesterday]
  people_o <- visits$participant_id[answered & !visits$yesterday]
  eligible <- answered & visits$participant_id %in% intersect(people_y, people_o)

  stages <- list(mbx_cd_visits = rep(TRUE, nrow(visits)),
                 mbx_positive_butyrate = positive,
                 mbx_valid_whole_grain = answered,
                 mbx_butyrate_both_states = eligible)
  flow <- data.frame(stage = names(stages), unit = "visit",
                     n = unname(vapply(stages, sum, integer(1))),
                     participants = unname(vapply(stages, function(keep) {
                       length(unique(visits$participant_id[keep]))
                     }, integer(1))))
  list(visits = visits[eligible, , drop = FALSE], flow = flow)
}

build_butyrate_cohort <- function(metadata_path, biom_path) {
  grain <- unname(food_fields[["diet_whole_grains"]])
  fields <- c("Participant ID", "External ID", "site_sub_coll", grain)
  cd <- read_study_metadata(metadata_path, fields, "metabolomics")
  visits <- data.frame(participant_id = cd[["Participant ID"]], sample_id = cd[["External ID"]],
                       collection_id = cd$site_sub_coll, whole_grain = cd[[grain]])

  # BIOM puts feature annotations after its sparse matrix. Use the location in
  # the checksum-pinned source to retain one row, then verify its ID and method.
  assay <- read_metabolomics(biom_path, visits$sample_id, feature_row = 41390L)
  feature <- assay$features
  peak <- "HILn_QI36"
  stopifnot(identical(feature$id, peak),
            identical(feature$metadata$Method, "HILIC-neg"),
            identical(feature$metadata$Metabolite, "butyrate"))
  butyrate_membership(visits, assay$intensity[, peak])
}
