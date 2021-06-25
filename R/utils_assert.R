assert_is <- function(x, what, name = deparse(substitute(x))) {
  if (!inherits(x, what)) {
    stop(sprintf("'%s' must be a %s", name,
                 paste(what, collapse = " / ")), call. = FALSE)
  }
  invisible(x)
}


assert_named <- function(x, unique = FALSE, name = deparse(substitute(x))) {
  if (is.null(names(x))) {
    stop(sprintf("'%s' must be named", name), call. = FALSE)
  }
  if (unique && any(duplicated(names(x)))) {
    stop(sprintf("'%s' must have unique names", name), call. = FALSE)
  }
}


assert_integer <- function(x, name = deparse(substitute(x)),
                           what = "integer") {
  if (!(is.integer(x))) {
    eps <- sqrt(.Machine$double.eps)
    usable_as_integer <- is.numeric(x) && (max(abs(round(x) - x)) < eps)
    if (!usable_as_integer) {
      stop(sprintf("'%s' must be an %s", name, what), call. = FALSE)
    }
    x <- as.integer(round(x))
  }
  invisible(x)
}

assert_positive_integer <- function(x, name = deparse(substitute(x))) {
    force(name)
    x <- assert_integer(x, name)
    if (any(x < 1L)) {
      stop(sprintf("'%s' must be at least 1", name), call. = FALSE)
    }
    invisible(x)
}


assert_character <- function(x, name = deparse(substitute(x))) {
  if (!(is.character(x))) {
    stop(sprintf("'%s' must be a character", name), call. = FALSE)
  }
  invisible(x)
}


assert_strictly_increasing <- function(x, name = deparse(substitute(x))) {
  if (any(diff(x) <= 0)) {
    stop(sprintf("'%s' must be strictly increasing", name), call. = FALSE)
  }
  invisible(x)
}


assert_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1L) {
    stop(sprintf("'%s' must be a scalar", name), call. = FALSE)
  }
  invisible(x)
}


assert_scalar_positive_integer <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  x <- assert_integer(x, name)
  if (x < 1L) {
    stop(sprintf("'%s' must be at least 1", name), call. = FALSE)
  }
  invisible(x)
}

assert_scalar_unit_interval <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  if (x < 0 | x > 1) {
    stop(sprintf("'%s' must be between 0 and 1", name), call. = FALSE)
  }
  invisible(x)
}

