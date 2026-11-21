food_fields <- c(
  diet_whole_grains = "Whole grains (wheat, oats, brown rice, rye, quinoa, wheat bread, wheat pasta)",
  diet_vegetables = "Vegetables (salad, tomatoes, onions, greens, carrots, peppers, green beans, etc)",
  diet_fruit = "Fruits (no juice) (Apples, raisins, bananas, oranges, strawberries, blueberries",
  diet_red_meat = "Red meat (beef, hamburger, pork, lamb)"
)

food_answers <- c(
  "No, I did not consume these products in the last 7 days",
  "Within the past 4 to 7 days", "Within the past 2 to 3 days",
  "Yesterday, 1 to 2 times", "Yesterday, 3 or more times"
)

build_reference_cohort <- function(path) {
  fields <- c("Participant ID", "site_sub_coll", "interval_days", "week_num",
              "diagnosis", "External ID", "data_type", "Antibiotics", unname(food_fields))
  headers <- names(read.csv(path, nrows = 0L, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  if (!all(fields %in% headers) || anyDuplicated(headers[headers %in% fields])) {
    stop("Missing or duplicate metadata field")
  }
  data <- read.csv(path, check.names = FALSE, na.strings = NULL, fileEncoding = "UTF-8-BOM",
                   colClasses = ifelse(headers %in% fields, "character", "NULL"))
  data <- data[data$diagnosis %in% c("CD", "nonIBD") &
                 data$data_type %in% c("metagenomics", "metabolomics"), ]
  keys <- c("Participant ID", "site_sub_coll", "interval_days")
  stopifnot(all(vapply(data[keys], function(x) all(nzchar(x)), logical(1))),
            all(is.finite(as.numeric(data$interval_days))))
  groups <- split(seq_len(nrow(data)), do.call(paste, c(data[keys], sep = "\t")))
  complete_foods <- Reduce(`&`, lapply(data[food_fields], function(x) x %in% food_answers))

  pairs <- lapply(groups, function(indices) {
    mgx <- indices[data$data_type[indices] == "metagenomics"]
    mbx <- indices[data$data_type[indices] == "metabolomics"]
    if (!length(mgx) || !length(mbx)) return(NULL)
    if (length(mbx) != 1L) stop("Multiple metabolomic profiles at a collection")
    group <- data[indices, ]
    ids <- data[["External ID"]]
    if (length(unique(group$diagnosis)) != 1L ||
        !all(group$Antibiotics %in% c("Yes", "No")) ||
        length(unique(group$Antibiotics)) != 1L ||
        !complete_foods[mbx] || !any(complete_foods[mgx])) return(NULL)

    agreement <- vapply(group[food_fields], function(x) {
      length(unique(x[x %in% food_answers])) == 1L
    }, logical(1))
    if (!all(agreement)) return(NULL)
    exact <- mgx[ids[mgx] == ids[mbx]]
    if (length(exact) == 1L) {
      microbe <- exact
      rule <- "exact_external_id"
    } else if (length(mgx) == 1L) {
      microbe <- mgx
      rule <- "unique_documentary_collection_key"
    } else {
      stop("Unresolved microbial profile choice at metabolomic ID ", ids[mbx])
    }
    week <- as.numeric(data$week_num[c(microbe, mbx)])
    stopifnot(complete_foods[microbe], all(is.finite(week)), week[1L] == week[2L])
    record <- data.frame(
      participant_id = data[["Participant ID"]][mbx],
      collection_id = data$site_sub_coll[mbx],
      interval_days = as.numeric(data$interval_days[mbx]),
      week_num = week[2L],
      diagnosis = data$diagnosis[mbx],
      mgx_sample_id = ids[microbe], mbx_sample_id = ids[mbx],
      pairing_rule = rule, mgx_candidate_count = length(mgx),
      antibiotics = data$Antibiotics[mbx]
    )
    record[names(food_fields)] <- unname(as.list(data[mbx, food_fields]))
    record
  })
  cohort <- do.call(rbind, unname(Filter(Negate(is.null), pairs)))
  cohort <- cohort[order(cohort$participant_id, cohort$interval_days, cohort$collection_id,
                         method = "radix"), ]
  rownames(cohort) <- NULL
  stopifnot(!anyDuplicated(cohort[c("participant_id", "collection_id", "interval_days")]),
            !anyDuplicated(cohort$mgx_sample_id), !anyDuplicated(cohort$mbx_sample_id),
            !anyDuplicated(unique(cohort[c("participant_id", "diagnosis")])$participant_id))
  cohort
}

read_genus <- function(path, ids) {
  source <- read.delim(path, check.names = FALSE, quote = "", comment.char = "")
  genus <- grepl("(^|\\|)g__[^|]+$", source[[1L]])
  names <- source[[1L]][genus]
  columns <- match(paste0(ids, "_profile"), colnames(source))
  stopifnot(length(names) == 187L, !anyDuplicated(names), !anyNA(columns))
  matrix <- t(as.matrix(source[genus, columns, drop = FALSE]))
  storage.mode(matrix) <- "double"
  dimnames(matrix) <- list(ids, names)
  stopifnot(all(is.finite(matrix)), all(matrix >= 0))
  matrix
}
