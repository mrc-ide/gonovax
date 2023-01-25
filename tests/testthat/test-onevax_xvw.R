context("onevax_xvw")

test_that("run_onevax_xvw works correctly", {
  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  y1 <- run_onevax_xvw(tt, gp, vea = 0, dur = 1e3)[[1]]

  # check no-one is offered vaccine with v switched off
  expect_true(all(y1$cum_vaccinated == 0))
  expect_true(all(y1$cum_offered == 0))

  y2 <- run_onevax_xvw(tt, gp, vea = 0, dur = 1e3, vbe = 1)[[1]]
  # check 100% vbe vaccinates all new entrants
  expect_equal(diff(rowSums(y2$cum_vaccinated[, , 1, ])), rep(12e3, max(tt)))
  # and no-one else
  expect_equal(sum(y2$cum_vaccinated[, , 2:3, ]), 0)
  # only adolescents are offered vaccine and all accept
  expect_equal(y2$cum_offered, y2$cum_vbe)

  # check can restart
  y2 <- run_onevax_xvw(tt, gp, vea = 0, dur = 1e3, vbe = 1)
  init_params <- lapply(y2, restart_params)
  y3 <- run_onevax_xvw(seq(max(tt), length.out = 2, by = 1),
                        gp, init_params,
                        vea = 0, dur = 1e3, vbe = 1)
  for (i in seq_along(y2)) {
    expect_equal(y2[[i]]$U[length(tt), , , ], y3[[i]]$U[1, , , ])
    expect_equal(y2[[i]]$I[length(tt), , , ], y3[[i]]$I[1, , , ])
    expect_equal(y2[[i]]$A[length(tt), , , ], y3[[i]]$A[1, , , ])
    expect_equal(y2[[i]]$S[length(tt), , , ], y3[[i]]$S[1, , , ])
    expect_equal(y2[[i]]$T[length(tt), , , ], y3[[i]]$T[1, , , ])
  }

  uptake <- c(0.5, 1) ## i.e. 50% for 1st param set, 100% for second
  # check VoD is working correctly
  y3e <- run_onevax_xvw(tt, gp, vea = 1, dur = 1e3, strategy = "VoD",
                         uptake = uptake)

  for (i in seq_along(y3e)) {
    # no-one in stratum V or W is vaccinated again
    expect_equal(sum(y3e[[i]]$cum_vaccinated[, , "V", ]), 0)
    expect_equal(sum(y3e[[i]]$cum_vaccinated[, , "W", ]), 0)
    # uptake % of treated are vaccinated
    expect_equal(y3e[[i]]$cum_vaccinated[, , "X", ] / uptake[i],
                 y3e[[i]]$cum_treated[, , "X", ])
    # all those treated are offered vaccine
    expect_equal(y3e[[i]]$cum_offered[, , "X", ],
                 y3e[[i]]$cum_treated[, , "X", ])
    # only those in X are offered vaccine
    expect_equal(sum(y3e[[i]]$cum_offered[, , c("V", "W"), ]), 0)
    # efficacy is perfect
    expect_equal(sum(y3e[[i]]$cum_treated[, , "V", ]), 0)
    # no-one is lost
    expect_equal(apply(y3e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }

  # check VoA is working correctly
  y4e <- run_onevax_xvw(tt, gp, vea = 1, dur = 1e3, strategy = "VoA",
                         uptake = uptake)
  for (i in seq_along(y4e)) {
  # no-one in stratum V or W is vaccinated again
  expect_equal(sum(y4e[[i]]$cum_vaccinated[, , "V", ]), 0)
  expect_equal(sum(y4e[[i]]$cum_vaccinated[, , "W", ]), 0)

  # uptake % of treated are vaccinated
  expect_equal(y4e[[i]]$cum_vaccinated[, , "X", ] / uptake[i],
             y4e[[i]]$cum_screened[, , "X", ] + y4e[[i]]$cum_treated[, , "X", ])

  # all those in X treated  or screened are offered vaccine
  expect_equal(y4e[[i]]$cum_offered[, , "X", ],
             y4e[[i]]$cum_treated[, , "X", ] + y4e[[i]]$cum_screened[, , "X", ])
  # only those in X are offered vaccine
  expect_equal(sum(y4e[[i]]$cum_offered[, , c("V", "W"), ]), 0)

  # efficacy is perfect
  expect_equal(sum(y4e[[i]]$cum_treated[, , "V", ]), 0)
  # no-one is lost
  expect_equal(apply(y4e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }

  # check vaccination targeting
  y5e <- run_onevax_xvw(tt, gp, vea = 1, dur = 1e3, strategy = "VoD(L)+VoA(H)",
                         uptake = uptake)
  for (i in seq_along(y5e)) {
    # no-one in stratum V or W is vaccinated again
    expect_equal(sum(y5e[[i]]$cum_vaccinated[, , "V", ]), 0)
    expect_equal(sum(y5e[[i]]$cum_vaccinated[, , "W", ]), 0)

    # only treated L are offered vaccine
    expect_equal(y5e[[i]]$cum_offered[, "L", "X", ],
                 y5e[[i]]$cum_treated[, "L", "X", ])
    # all attending H are offered vaccine
    expect_equal(y5e[[i]]$cum_offered[, "H", "X", ],
                 y5e[[i]]$cum_treated[, "H", "X", ] +
                   y5e[[i]]$cum_screened[, "H", "X", ])

    # check that vaccinated = offered * uptake
    expect_equal(y5e[[i]]$cum_offered * uptake[i], y5e[[i]]$cum_vaccinated)

    # only those in X are offered vaccine
    expect_equal(sum(y4e[[i]]$cum_offered[, , c("V", "W"), ]), 0)

    # efficacy is perfect
    expect_equal(sum(y4e[[i]]$cum_treated[, , "V", ]), 0)
    # no-one is lost
    expect_equal(apply(y4e[[i]]$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)
  }

  # check length of uptake vector must be 1 or length(gp)
  expect_error(run_onevax_xvw(tt, gp, vea = 1, dur = 1e3, strategy = "VbE",
                  uptake = c(0, 0.5, 1)))
  expect_error(run_onevax_xvw(tt, gp, vea = c(0, 1, 2), dur = 1e3,
                               strategy = "VbE", uptake = 1))
  expect_error(run_onevax_xvw(tt, gp, vea = 1, dur = c(0, 1e2, 1e3),
                               strategy = "VbE", uptake = 1))
})

test_that("vaccine effects work as expected", {
  tt <- seq(0, 5)
  gp <- gono_params(1:2)

  ## check perfect protection against symptoms works
  y1 <- run_onevax_xvw(tt, gp, ves = 1, dur = 1,
                        uptake = 1, strategy = "VoA")[[1]]    # if everyone takes the vaccine, and it protects from symptoms...

  expect_true(all(y1$S[, , 2, ] == 0))           #vaccinated people with symptoms are zero
  expect_true(all(y1$S[-1, , 3, ] > 1e-6))       #waned people with symptoms are greater than zero because they lose protection
  expect_true(all(y1$A[-1, , , ] > 1e-6))        #the total thos asymptomatic is greater than zero across time and categories 

  ## check perfect protection against duration works
  y2 <- run_onevax_xvw(tt, gp, ved = 1, dur = 1,
                        uptake = 1, strategy = "VoA")[[1]]

  expect_true(all(y2$A[, , 2, ] < 1e-6))
  expect_true(all(y2$A[-1, , 3, ] > 1e-6))
  expect_true(all(y2$S[-1, , , ] > 0))

})

test_that("can set initial coverage", {
  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  cov <- 0.33

  ## check perfect protection against symptoms works
  y <- run_onevax_xvw(tt, gp, ves = 1, dur = 1e3, coverage = 0.33)[[1]]
  expect_equivalent(apply(y$N[, , -1], c(1, 2), sum) / apply(y$N, c(1, 2), sum),
                    rep(cov, length(tt) * 2))

  expect_error(run_onevax_xvw(tt, gp, ves = 1, dur = 1e3, coverage = c(0.1, 0)),
               "'coverage' must be a scalar")
})
