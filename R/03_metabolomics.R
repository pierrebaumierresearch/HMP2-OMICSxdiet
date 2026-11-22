# Read the hash-pinned sparse BIOM 1.0 source in chunks. Parsing the full JSON
# array at once would create millions of temporary R objects.
read_metabolomics <- function(path, ids, feature_row = NULL) {
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  chunk_size <- 4L * 1024L^2L
  read_chunk <- function() rawToChar(readBin(con, "raw", n = chunk_size))
  first <- read_chunk()
  marker <- '"data": ['
  start <- regexpr(marker, first, fixed = TRUE)[1L]
  stopifnot(start > 0L)
  prefix <- paste0(substr(first, 1L, start - 1L), '"data": []')
  header <- jsonlite::fromJSON(paste0(prefix, "}"), simplifyVector = FALSE)
  stopifnot(identical(header$matrix_type, "sparse"),
            identical(unlist(header$shape), c(81867L, 546L)))
  selected_rows <- seq_len(81867L)
  if (!is.null(feature_row)) {
    stopifnot(is.numeric(feature_row), length(feature_row) == 1L, feature_row %in% selected_rows)
    selected_rows <- feature_row
  }
  # Linear index matches source order, with samples along the first dimension.
  values <- matrix(0, nrow = 546L, ncol = length(selected_rows))
  buffer <- substr(first, start + nchar(marker), nchar(first))
  last_position <- 0
  entries <- 0
  repeat {
    end <- regexpr("]]", buffer, fixed = TRUE)[1L]
    if (end > 0L) {
      boundary <- end
    } else {
      separators <- gregexpr("],", buffer, fixed = TRUE)[[1L]]
      boundary <- tail(separators, 1L)
      if (boundary < 1L) stop("Truncated or malformed BIOM data")
    }
    tokens <- chartr("[],", "   ", substr(buffer, 1L, boundary))
    triplets <- scan(text = tokens, what = double(), quiet = TRUE)
    stopifnot(length(triplets) > 0L, length(triplets) %% 3L == 0L)
    i <- seq.int(1L, length(triplets), by = 3L)
    feature_index <- triplets[i]
    sample_index <- triplets[i + 1L]
    intensity <- triplets[i + 2L]
    stopifnot(all(feature_index == trunc(feature_index) & feature_index >= 0 & feature_index < 81867L),
              all(sample_index == trunc(sample_index) & sample_index >= 0 & sample_index < 546L),
              all(is.finite(intensity) & intensity >= 0))
    positions <- feature_index * 546L + sample_index + 1L
    stopifnot(positions[1L] > last_position, !is.unsorted(positions, strictly = TRUE))
    if (is.null(feature_row)) {
      values[positions] <- intensity
    } else {
      keep <- feature_index == feature_row - 1L
      values[sample_index[keep] + 1L, 1L] <- intensity[keep]
    }
    last_position <- tail(positions, 1L)
    entries <- entries + length(positions)
    if (end > 0L) break
    next_chunk <- read_chunk()
    if (!nzchar(next_chunk)) stop("Missing BIOM data terminator")
    buffer <- paste0(substr(buffer, boundary + 2L, nchar(buffer)), next_chunk)
  }

  suffix <- substr(buffer, end + 2L, nchar(buffer))
  repeat {
    part <- read_chunk()
    if (!nzchar(part)) break
    suffix <- paste0(suffix, part)
  }
  metadata <- jsonlite::fromJSON(paste0(prefix, suffix))
  stopifnot(nrow(metadata$columns) == 546L, nrow(metadata$rows) == 81867L,
            !anyDuplicated(metadata$columns$id), !anyDuplicated(metadata$rows$id))
  selected <- match(ids, metadata$columns$id)
  stopifnot(!anyNA(selected), entries == 27348904)
  matrix <- values[selected, , drop = FALSE]
  features <- metadata$rows[selected_rows, , drop = FALSE]
  dimnames(matrix) <- list(ids, features$id)
  list(intensity = matrix, features = features)
}
