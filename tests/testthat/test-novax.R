test_that("run_novax works correctly", {
  tt <- c(0, 1)
  y1 <- run_novax_int(1, tt)
  y2 <- run_novax_int(2, tt)

  y1e <- run_novax_int(1, tt, equilib = TRUE)
  y2e <- run_novax_int(2, tt, equilib = TRUE)

  expect_error(run_novax_int(1:2, tt, equilib = TRUE),
               label = "if length(n) > 1, equilib must be FALSE")

  y1m <- run_novax(1, tt)
  y12m <- run_novax(1:2, tt)
  y1me <- run_novax(1, tt, equilib = TRUE)
  y12me <- run_novax(1:2, tt, equilib = TRUE)

  expect_equal(y1, y1m[[1]])
  expect_equal(y1e, y1me[[1]])
  expect_equal(y1, y12m[[1]])
  expect_equal(y2, y12m[[2]])
  expect_equal(y1e, y12me[[1]])
  expect_equal(y2e, y12me[[2]])

})
