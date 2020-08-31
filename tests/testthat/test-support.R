testthat("aggregate works", {
  dt <- 1 / 12
  tt <- seq(0, 1, dt)
  y <- run_onevax(1:2, tt, eff = 1, dur = 4)
  z <- aggregate(y, "cum_incid", FALSE)
  expect_equal(dim(z), c(2, length(tt)))
  expect_equal(z[1, ], apply(y[[1]]$cum_incid, 1, sum))
  expect_equal(z[2, ], apply(y[[2]]$cum_incid, 1, sum))

  z1 <- aggregate(y, "cum_incid", TRUE)
  expect_equal(z1, t(apply(z, 1, diff)) / dt)

  z2 <- aggregate(y, "cum_incid", FALSE, mean)
  expect_equal(z2, colMeans(z))

  z3 <- aggregate(y, "cum_incid", TRUE, mean)
  expect_equal(z3, colMeans(z1))

})