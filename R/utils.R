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

switch_levels <- function(x) {
  nms <- names(x[[1]])
  y <- lapply(nms, function(z) lapply(x, "[[", z))
  names(y) <- nms
  y
}

vlapply <- function(x, fun, ...) {
  vapply(x, fun, logical(1), ...)
}

vnapply <- function(x, fun, ...) {
  vapply(x, fun, numeric(1), ...)
}
