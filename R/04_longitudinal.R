# Select assay-specific Crohn's records using exact exported field names.
read_study_metadata <- function(path, fields, assay_type) {
  fields <- unique(c("diagnosis", "data_type", fields))
  headers <- names(read.csv(path, nrows = 0L, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  if (!all(fields %in% headers) || anyDuplicated(headers[headers %in% fields])) {
    stop("Missing or duplicate study metadata field")
  }
  records <- read.csv(path, check.names = FALSE, na.strings = NULL, fileEncoding = "UTF-8-BOM",
                      colClasses = ifelse(headers %in% fields, "character", "NULL"))
  records[records$diagnosis == "CD" & records$data_type == assay_type, , drop = FALSE]
}

# Indices refer to the input order. Form neighbors before checking either endpoint
# so that an ambiguous day, an unusable profile or an unanswered kit breaks the pair.
adjacent_kit_rows <- function(participant_id, receipt_date, week_num, eligible) {
  n <- length(participant_id)
  stopifnot(length(receipt_date) == n, length(week_num) == n, length(eligible) == n,
            inherits(receipt_date, "Date"), !anyNA(receipt_date), !anyNA(participant_id))
  day <- data.frame(participant_id, receipt_date)
  single <- !duplicated(day) & !duplicated(day, fromLast = TRUE)
  idx <- order(participant_id, receipt_date, method = "radix")
  a <- head(idx, -1L)
  b <- tail(idx, -1L)
  gap <- as.integer(receipt_date[b] - receipt_date[a])
  keep <- single[a] & single[b] & eligible[a] & eligible[b] &
    participant_id[a] == participant_id[b] & gap >= 7L & gap <= 28L &
    is.finite(week_num[a]) & is.finite(week_num[b]) & week_num[b] > week_num[a]
  data.frame(a = a[which(keep)], b = b[which(keep)])
}

# Both a switch and a stable interval must exist for this person and starting state.
same_start_support <- function(participant_id, start, end, eligible) {
  stopifnot(is.logical(start), is.logical(end), is.logical(eligible),
            length(participant_id) == length(start), length(end) == length(start),
            length(eligible) == length(start), !anyNA(participant_id))
  usable <- eligible %in% TRUE & !is.na(start) & !is.na(end)
  group <- interaction(participant_id, start, drop = TRUE)
  supported <- intersect(group[usable & start != end], group[usable & start == end])
  usable & group %in% supported
}

# Assign on the complete set of source participants before joining eligible intervals.
participant_folds <- function(ids) {
  stopifnot(length(ids) > 0L, !anyNA(ids), all(nzchar(ids)))
  old_locale <- Sys.getlocale("LC_COLLATE")
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    Sys.setlocale("LC_COLLATE", old_locale)
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  })
  if (!nzchar(Sys.setlocale("LC_COLLATE", "C.UTF-8"))) stop("C.UTF-8 collation is required")
  people <- sort(unique(ids))
  set.seed(20261121L, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  folds <- sample(rep(seq_len(4L), length.out = length(people)))
  folds[match(ids, people)]
}

build_mgx_cohorts <- function(metadata_path, taxon_path) {
  dictionary <- read.csv("data/diet_items.csv", colClasses = "character")
  stopifnot(nrow(dictionary) == 21L, !anyDuplicated(dictionary$item),
            !anyDuplicated(dictionary$source_field), "whole_grains" %in% dictionary$item)
  med_screen <- "1) In the past 2 weeks, have you received any of the following medications"
  fields <- c("Participant ID", "External ID", "week_num", "date_of_receipt",
              "Antibiotics", med_screen, dictionary$source_field)
  m <- read_study_metadata(metadata_path, fields, "metagenomics")
  person <- m[["Participant ID"]]
  id <- m[["External ID"]]
  receipt <- as.Date(m$date_of_receipt, format = "%Y-%m-%d")
  week <- as.numeric(m$week_num)
  stopifnot(length(id) > 0L, all(nzchar(id)), all(nzchar(person)), !anyDuplicated(id),
            !anyNA(receipt), all(format(receipt, "%Y-%m-%d") == m$date_of_receipt))
  day <- data.frame(person, receipt)
  single <- !duplicated(day) & !duplicated(day, fromLast = TRUE)

  taxon_headers <- names(read.delim(taxon_path, nrows = 0L, check.names = FALSE,
                                    quote = "", comment.char = ""))
  stopifnot(!anyDuplicated(taxon_headers))
  has_profile <- paste0(id, "_profile") %in% taxon_headers
  abundance <- read_genus(taxon_path, id[has_profile])
  positive <- has_profile
  positive[has_profile] <- rowSums(abundance) > 0
  valid_day <- single & positive
  if (!any(valid_day)) stop("No unambiguous day with a positive microbial profile")
  # The prevalence denominator includes valid days whose dietary recall is missing.
  detected <- colSums(abundance[id[valid_day], , drop = FALSE] > 0)
  panel <- names(detected)[2L * detected >= sum(valid_day)]
  if (!length(panel)) stop("No genus meets the common-panel prevalence threshold")
  panel <- panel[order(sub("^.*[|]g__", "", panel), method = "radix")]
  panel_day <- has_profile
  panel_day[has_profile] <- rowSums(abundance[, panel, drop = FALSE]) > 0

  diet <- m[dictionary$source_field]
  names(diet) <- dictionary$item
  grain <- diet$whole_grains
  valid_recall <- valid_day & grain %in% food_answers
  all_items <- Reduce(`&`, lapply(diet, function(x) x %in% food_answers))
  yesterday <- grain %in% c("Yesterday, 1 to 2 times", "Yesterday, 3 or more times")
  three_days <- yesterday | grain == "Within the past 2 to 3 days"
  known_negative <- m[[med_screen]] %in% c("Yes", "No") & m$Antibiotics %in% "No"
  fold <- participant_folds(person)

  pair <- adjacent_kit_rows(person, receipt, week, valid_recall)
  a <- pair$a
  b <- pair$b
  intervals <- data.frame(
    participant_id = person[a], earlier_id = id[a], later_id = id[b],
    gap_days = as.integer(receipt[b] - receipt[a]), outer_fold = fold[a],
    start_y = yesterday[a], end_y = yesterday[b],
    start_three_day = three_days[a], end_three_day = three_days[b],
    start_screen = m[[med_screen]][a], end_screen = m[[med_screen]][b],
    start_antibiotic = m$Antibiotics[a], end_antibiotic = m$Antibiotics[b],
    both_antibiotic_negative = known_negative[a] & known_negative[b],
    panel_eligible = panel_day[a] & panel_day[b], whole_recall = all_items[a] & all_items[b]
  )
  all_pairs <- rep(TRUE, nrow(intervals))
  support <- function(start, end, eligible) same_start_support(person[a], start, end, eligible)
  full_support <- support(intervals$start_y, intervals$end_y, all_pairs)
  intervals$primary <- support(intervals$start_y, intervals$end_y, intervals$panel_eligible)
  intervals$three_day <- support(intervals$start_three_day, intervals$end_three_day,
                                 intervals$panel_eligible)
  intervals$antibiotic_negative <- support(intervals$start_y, intervals$end_y,
                                            intervals$panel_eligible & intervals$both_antibiotic_negative)

  flow_row <- function(stage, unit, ids) {
    data.frame(stage = paste0("mgx_", stage), unit = unit, n = length(ids),
               participants = length(unique(ids)))
  }
  observed <- !duplicated(day)
  flow <- rbind(flow_row("source", "assay", person), flow_row("observed_days", "day", person[observed]),
                flow_row("valid_profile", "day", person[valid_day]),
                flow_row("valid_recall", "day", person[valid_recall]))
  stages <- list(all_intervals = all_pairs, same_start_before_panel = full_support,
                 panel_intervals = intervals$panel_eligible, primary = intervals$primary,
                 primary_toward_switch = intervals$primary & !intervals$start_y & intervals$end_y,
                 primary_toward_stable = intervals$primary & !intervals$start_y & !intervals$end_y,
                 primary_away_switch = intervals$primary & intervals$start_y & !intervals$end_y,
                 primary_away_stable = intervals$primary & intervals$start_y & intervals$end_y,
                 three_day = intervals$three_day, antibiotic_negative = intervals$antibiotic_negative,
                 whole_recall = intervals$whole_recall,
                 whole_recall_antibiotic_negative = intervals$whole_recall & intervals$both_antibiotic_negative)
  for (stage in names(stages)) {
    flow <- rbind(flow, flow_row(stage, "interval", intervals$participant_id[stages[[stage]]]))
  }
  list(intervals = intervals, panel = panel, flow = flow)
}
