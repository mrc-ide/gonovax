test_that("aggregate works", {
  dt <- 1 / 12
  tt <- seq(0, 1, dt)
  y <- run_onevax(1:2, tt, eff = 1, dur = 4, vd = 1)
  z <- aggregate(y, "cum_incid", FALSE)
  expect_equal(dim(z), c(2, length(tt)))
  expect_equal(z[1, ], apply(y[[1]]$cum_incid, 1, sum))
  expect_equal(z[2, ], apply(y[[2]]$cum_incid, 1, sum))

  z1 <- aggregate(y, "cum_incid", TRUE)
  expect_equal(z1, t(apply(z, 1, diff)) / dt)

  z2 <- aggregate(y, "cum_incid", FALSE, f = mean)
  expect_equal(z2, colMeans(z))

  z3 <- aggregate(y, "cum_incid", TRUE, f = mean)
  expect_equal(z3, colMeans(z1))

  z4 <- aggregate(y, "cum_incid", TRUE, stratum = 1, f = mean)
  z5 <- aggregate(y, "cum_incid", TRUE, stratum = 3, f = mean)
  expect_equal(z4 + z5, z3)

})

test_that("extract_flows works", {
  tt <- seq(0, 2)
  y <- run_onevax(1:2, tt, eff = 1, dur = 4, vd = 1, ve = 1)
  z <- extract_flows(y)

  expect_equal(z$cum_incid[1, ], z$incid[1, ])
  expect_equal(z$cum_incid[2, ] - z$cum_incid[1, ], z$incid[2, ])
  expect_equal(z$vaccinated + z$revaccinated,
               t(aggregate(y, "cum_vaccinated", as_incid = TRUE)))
  z$vaccinated
  z$revaccinated
})
