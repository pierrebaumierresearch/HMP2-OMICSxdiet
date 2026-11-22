sha256 <- function(path) {
  command <- Sys.which("sha256sum")
  if (!nzchar(command)) stop("sha256sum is required. Use WSL on Windows.")
  output <- system2(command, c("--", shQuote(path)), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(output, "status")) || length(output) != 1L ||
      !grepl("^[0-9a-f]{64} ", output)) stop("Unable to hash ", path)
  substr(output, 1L, 64L)
}

verified_sources <- function(directory, download = FALSE) {
  manifest <- read.csv("data/sources.csv", colClasses = "character")
  if (download) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  paths <- setNames(file.path(directory, manifest$file), manifest$file)
  for (i in seq_along(paths)) {
    path <- paths[[i]]
    if (!file.exists(path) && download) {
      pending <- tempfile(".download-", tmpdir = directory)
      tryCatch({
        status <- download.file(manifest$url[i], pending, mode = "wb")
        if (status != 0L || !identical(sha256(pending), manifest$sha256[i])) {
          stop("Source download failed its checksum: ", path)
        }
        if (!file.rename(pending, path)) stop("Could not save ", path)
      }, finally = unlink(pending))
    }
    if (!file.exists(path)) stop("Missing input: ", path, ". Run with --download to acquire it.")
    if (!identical(sha256(path), manifest$sha256[i])) stop("Source checksum mismatch: ", path)
  }
  paths
}

study_inputs <- function(entry_point, args = commandArgs(trailingOnly = TRUE)) {
  if (!identical(args, character()) && !identical(args, "--download")) {
    stop("Usage from the project root: Rscript ", entry_point, " [--download]")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Missing jsonlite. Restore the project environment first.")
  }
  options(timeout = max(300, getOption("timeout")))
  raw_dir <- Sys.getenv("HMP2_RAW_DIR", unset = "data/raw")
  verified_sources(raw_dir, download = identical(args, "--download"))
}

write_output <- function(path, write_file) {
  parent <- dirname(path)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  pending <- tempfile(".output-", tmpdir = parent)
  on.exit(unlink(pending), add = TRUE)
  write_file(pending)
  if (!file.rename(pending, path)) stop("Could not save ", path)
}
