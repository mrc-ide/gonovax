`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}

gonovax_file <- function(...) {
  system.file(..., package = "gonovax", mustWork = TRUE)
}

read_csv <- function(...) {
  utils::read.csv(..., check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv <- function(...) {
  utils::write.csv(..., row.names = FALSE)
}

data_frame <- function(...) {
  data.frame(..., check.names = FALSE, stringsAsFactors = FALSE)
}

