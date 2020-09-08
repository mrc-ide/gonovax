context("utils")

test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})

test_that("data_frame and read_csv do not do factor conversion", {
  d <- data_frame(a = letters, b = LETTERS, c = 1:26)
  expect_type(d$a, "character")
  expect_type(d$b, "character")
  path <- tempfile()
  on.exit(unlink(path))
  write_csv(d, path)
  expect_equal(read_csv(path), d)
})

test_that("sircovid_file throws for missing files", {
  expect_true(file.exists(gonovax_file("odin/model.R")))
  ## NOTE: not testing error string because it comes from base R
  expect_error(sircovid_file("odin/acidic.R"))
})

test_that("switch_levels works as expected", {
  x <- list(a = list(b = 1, c = 2), d = list(b = 3, c = 4),
            e = list(b = 5, c = 6))
  y <- switch_levels(x)
  z <- switch_levels(y)
  expect_equal(y, list(b = list(a = 1, d = 3, e = 5),
                       c = list(a = 2, d = 4, e = 6)))
  expect_equal(x, z)
})
