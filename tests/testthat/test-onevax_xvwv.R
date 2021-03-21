context("onevax_xvwv")

test_that("run_onevax works correctly", {
  tt <- seq(0, 5)
  gp <- gono_params(1:2)
  y1 <- run_onevax_xvwv(tt, gp, eff = 0, dur = 1e3)[[1]]

  # check no-one is vaccinated with v switched off
  expect_true(all(y1$cum_vaccinated == 0))
  y2 <- run_onevax_xvwv(tt, gp, eff = 0, dur = 1e3, ve = 1)
  # check 100% vbe vaccinates all new entrants
  expect_equal(diff(rowSums(y2[[1]]$cum_vaccinated[, , 1])), rep(12e3, max(tt)))
  # and no-one else
  expect_equal(sum(y2[[1]]$cum_vaccinated[, , 2:3]), 0)

  # check can restart
  init_params <- lapply(y2, restart_params)
  y3 <- run_onevax_xvwv(seq(max(tt), length.out = 2, by = 1),
                        gp, init_params,
                        eff = 0, dur = 1e3, ve = 1)
  for (i in seq_along(y2)) {
    expect_equal(y2[[i]]$U[length(tt), , ], y3[[i]]$U[1, , ])
    expect_equal(y2[[i]]$I[length(tt), , ], y3[[i]]$I[1, , ])
    expect_equal(y2[[i]]$A[length(tt), , ], y3[[i]]$A[1, , ])
    expect_equal(y2[[i]]$S[length(tt), , ], y3[[i]]$S[1, , ])
    expect_equal(y2[[i]]$T[length(tt), , ], y3[[i]]$T[1, , ])
  }


  # check vaccination on treatment is working correctly
  y3e <- run_onevax_xvwv(tt, gp, eff = 1, dur = 1e3, vd = 1)[[1]]
  expect_equal(y3e$cum_vaccinated[, , 1], y3e$cum_treated[, , 1])
  expect_equal(y3e$cum_vaccinated[, , 3], y3e$cum_treated[, , 3])
  expect_equal(sum(y3e$cum_treated[, , 2]), 0)
  expect_equal(apply(y3e$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)

  # check vaccination on screening is working correctly
  y4e <- run_onevax_xvwv(tt, gp, eff = 1, dur = 1e3, vs = 1)[[1]]
  expect_equal(y4e$cum_vaccinated[, , 1], y4e$cum_screened[, , 1])
  expect_equal(sum(y4e$cum_treated[, , 2]), 0)
  expect_equal(y4e$cum_vaccinated[, , 3], y4e$cum_screened[, , 3])
  expect_equal(apply(y4e$N, 1, sum), rep(6e5, 6), tolerance = 1e-5)

  # check vaccination targeting
  v <- list(c(0, 1))
  y5e <- run_onevax_xvwv(tt, gp, eff = 1, dur = 1e3, vs = v, vd = v)
  for (i in seq_along(y5e)) {
    expect_equal(y5e[[i]]$N[, 1, 2], rep(0, 6))
    expect_equal(sum(y5e$cum_vaccinated[, 1, ]), 0)
    expect_true(all(y5e[[i]]$cum_vaccinated[-1, 2, 1] > 0))
    expect_true(all(y5e[[i]]$cum_vaccinated[-1, 2, 3] > 0))
    expect_true(all(y5e[[i]]$cum_vaccinated[1, 2, 2] == 0))
    expect_true(all(y5e[[i]]$N[-1, 2, 2] > 0))
  }
})
