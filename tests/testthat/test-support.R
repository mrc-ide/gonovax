test_that("aggregate works", {
  dt <- 1 / 12
  tt <- seq(0, 1, dt)
  y <- run_onevax_xvwv(tt, gono_params(1:2), eff = 1, dur = 4, uptake = 1,
                       strategy = "VoD", ve = 1)
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
  y <- run_onevax_xvwv(tt, gono_params(1:2), eff = 1, dur = 4, uptake = 1,
                       strategy = "VoD", ve = 1)
  z <- extract_flows(y)

  expect_equal(z$cum_treated[1, ], z$treated[1, ])
  expect_equal(z$cum_treated[2, ] - z$cum_treated[1, ], z$treated[2, ])
  expect_equal(z$vaccinated, t(aggregate(y, "cum_vaccinated", as_incid = TRUE)))
  expect_equal(z$revaccinated[2, ],
               sapply(y, function(x) diff(rowSums(x$cum_vaccinated[-1, , 3]))))
})

test_that("gonovax_year works as expected", {
  expect_equal(gonovax_year(2009), 0)
  expect_error(gonovax_year(2006),
               "Negative dates, gonovax_year likely applied twice")
  expect_equal(gonovax_year_as_year(2), 2011)
  expect_error(gonovax_year_as_year(-1), "'gonovax_year' must be at least 1")
})


test_that("dbetabinom", {

  ## when prob = 1/2, rho = 1/3 (equivalently a = b = 1),
  ## equivalent to discrete uniform from 0 to size
  expect_equal(dbetabinom(4, 35, 1 / 2, 1 / 3), 1 / 36)
  expect_equal(dbetabinom(4, 35, 1 / 2, 1 / 3, log = TRUE), log(1 / 36))


  ## compare with extraDistr::dbbinom (uses a and b as parameters - to match
  ## use prob = a / (a + b) and rho = 1 / (a + b + 1))
  f <- function(x, size, a, b, log = FALSE) {
    prob  <- a / (a + b)
    rho <- 1 / (a + b + 1)

    dbetabinom(x, size, prob, rho, log)
  }

  ## compare to extraDistr::dbbinom(15, 63, 2, 5) ~ 0.03613356
  expect_equal(f(15, 63, 2, 5), 0.03613356, tolerance = 1e-8)
  ## compare to extraDistr::dbbinom(15, 63, 2, 5, log = TRUE) ~ -3.320533
  expect_equal(f(15, 63, 2, 5, log = TRUE), -3.320533, tolerance = 1e-6)

  ## compare to extraDistr::dbbinom(672, 50454, 3, 2) ~ 4.18089e-08
  expect_equal(f(672, 50454, 3, 2), 4.18089e-08, tolerance = 1e-13)
  ## compare to extraDistr::dbbinom(15, 63, 2, 5, log = TRUE) ~ -16.99016
  expect_equal(f(672, 50454, 3, 2, log = TRUE), -16.99016, tolerance = 1e-5)
})
