# Run from the repository root with Rscript --vanilla tests/test_study_cohorts.R.
source("R/02_cohort.R")
source("R/04_longitudinal.R")
source("R/05_butyrate.R")

test_metagenomic_cohorts <- function() {
  fields <- read.csv("data/diet_items.csv", stringsAsFactors = FALSE)
  tea <- "Tea or coffee no sugar and no sugar replacement"
  grain <- unname(food_fields[["diet_whole_grains"]])
  screen <- "1) In the past 2 weeks, have you received any of the following medications"
  stopifnot(nrow(fields) == 21L, !anyDuplicated(fields$source_field),
            all(c(tea, grain) %in% fields$source_field),
            !any(c("Antibiotics", screen) %in% fields$source_field))
  kits <- function(person, days, responses) {
    id <- paste0(person, seq_along(days))
    rows <- data.frame(
      "Participant ID" = person, "External ID" = id, site_sub_coll = id,
      date_of_receipt = as.Date("2020-01-01") + days, week_num = days / 7,
      interval_days = days, diagnosis = "CD", data_type = "metagenomics",
      Antibiotics = "No", check.names = FALSE
    )
    rows[[screen]] <- "No"
    for (field in fields$source_field) rows[[field]] <- food_answers[1L]
    rows[[grain]] <- responses
    rows
  }
  o <- food_answers[3L]
  y <- food_answers[4L]
  visits <- rbind(
    kits("A", c(0, 7, 14, 21, 28), c(o, o, y, food_answers[5L], o)),
    kits("B", c(0, 7), c(o, o)),
    kits("C", c(0, 7), c(o, y)),
    kits("D", c(0, 7, 7, 14, 21, 49, 55, 84, 91), c(o, o, o, y, y, y, y, y, y)),
    kits("E", c(0, 7, 14, 21, 28), c(o, "", o, y, o)),
    kits("F", c(0, 7, 14, 21), c(o, o, y, o))
  )
  ids <- visits[["External ID"]]
  visits[[screen]][ids == "A1"] <- ""
  visits[[tea]][ids == "A2"] <- ""
  visits$week_num[ids == "D9"] <- visits$week_num[ids == "D8"]
  profiles <- matrix(0, nrow = 187L, ncol = length(ids),
                     dimnames = list(NULL, paste0(ids, "_profile")))
  profiles[1:22, ] <- 1
  profiles[, ids %in% c("E4", "F1")] <- 0
  profiles[23L, ids == "F1"] <- 1
  taxonomy <- paste0("k__Bacteria|g__Genus", seq_len(nrow(profiles)))
  files <- tempfile(c("metadata-", "taxa-"))
  on.exit(unlink(files))
  # Reverse the metadata order to exercise both the time sort and the assay-ID join.
  write.csv(visits[nrow(visits):1L, ], files[1L], row.names = FALSE)
  write.table(data.frame(taxon = taxonomy, profiles, check.names = FALSE),
              files[2L], sep = "\t", quote = FALSE, row.names = FALSE)
  result <- build_mgx_cohorts(files[1L], files[2L])
  pairs <- result$intervals
  keys <- paste(pairs$earlier_id, pairs$later_id, sep = "-")
  a <- c("A1-A2", "A2-A3", "A3-A4", "A4-A5")
  # D has a duplicate day, an eligible 28-day gap, gaps of 6 and 29 days,
  # and a final gap of 7 days without an increase in nominal study week.
  expected <- c(a, "B1-B2", "C1-C2", "D4-D5", "D5-D6", "F1-F2", "F2-F3", "F3-F4")
  stopifnot(!anyDuplicated(keys), setequal(keys, expected),
            setequal(result$panel, taxonomy[1:22]))
  # B and C supply opposite arms in different people. F loses the stable
  # comparator at the panel filter. Only A retains both arms for each start.
  stopifnot(setequal(keys[pairs$primary], a),
            setequal(keys[!pairs$panel_eligible], "F1-F2"))
  # Missing tea invalidates the wider recall only. A blank medication parent
  # makes the exported "No" unusable, removing A's only O-start stable arm.
  stopifnot(setequal(keys[!pairs$whole_recall], a[1:2]),
            setequal(keys[pairs$antibiotic_negative], a[3:4]))
  # Folds are assigned to all six source people, including E without intervals.
  expected_fold <- c(A = 4L, B = 2L, C = 1L, D = 2L, E = 3L, F = 1L)
  stopifnot(all(pairs$outer_fold == unname(expected_fold[pairs$participant_id])))
}

test_butyrate_cohort <- function() {
  visits <- data.frame(
    participant_id = c("A", "A", "B", "B", "C", "C"),
    sample_id = paste0("s", 1:6), collection_id = paste0("c", 1:6),
    whole_grain = food_answers[c(4L, 2L, 4L, 2L, 4L, 1L)]
  )
  # B and C each lose their O visit for a different reason.
  intensities <- setNames(c(2, 1, 3, 0, 4, NA_real_), visits$sample_id)
  result <- butyrate_membership(visits, rev(intensities))
  stopifnot(setequal(result$visits$sample_id, c("s1", "s2")),
            identical(result$visits$yesterday, c(TRUE, FALSE)))
  duplicate <- rbind(visits, visits[1L, ])
  stopifnot(inherits(try(butyrate_membership(duplicate, intensities), silent = TRUE), "try-error"),
            inherits(try(butyrate_membership(visits, intensities[-1L]), silent = TRUE), "try-error"))
}

test_metagenomic_cohorts()
test_butyrate_cohort()
message("Synthetic cohort checks passed.")
