context("onevax")

test_that("run_onevax works correctly", {
  tt <- seq(0, 5)
  gp <- gono_params(2)
  y1 <- run_onevax_xvwv(tt, gp, eff = 0, dur = 1e3)[[1]]

  # check no-one is vaccinated with v switched off
  expect_true(all(y1$cum_vaccinated == 0))
  y2 <- run_onevax_xvwv(tt, gp, eff = 0, dur = 1e3, ve = 1)[[1]]
  # check 100% vbe vaccinates all new entrants
  expect_equal(diff(rowSums(y2$cum_vaccinated[, , 1])), rep(12e3, max(tt)))
  # and no-one else
  expect_equal(sum(y2$cum_vaccinated[, , 2:3]), 0)

  # check can restart
  init_params <- restart_params(y2)
  y3 <- run_onevax_xvwv(seq(max(tt), length.out = 2, by = 1),
                        gp, list(init_params),
                        eff = 0, dur = 1e3, ve = 1)[[1]]
  plot(y2$t, apply(y2$U, 1, sum), type = "l", xlim = c(0, 10))
  lines(y3$t, apply(y3$U, 1, sum), col = "red")
  expect_equal(y2$U[length(y2$t), , ], y3$U[1, , ])
  expect_equal(y2$I[length(y2$t), , ], y3$I[1, , ])
  expect_equal(y2$A[length(y2$t), , ], y3$A[1, , ])
  expect_equal(y2$S[length(y2$t), , ], y3$S[1, , ])
  expect_equal(y2$T[length(y2$t), , ], y3$T[1, , ])

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
  y5e <- run_onevax_xvwv(tt, gp, eff = 1, dur = 1e3,
                    vs = c(0, 1), vd = c(0, 1))[[1]]
  expect_equal(y5e$N[, 1, 2], rep(0, 6))
  expect_equal(sum(y5e$cum_vaccinated[, 1, ]), 0)
  expect_true(all(y5e$cum_vaccinated[-1, 2, 1] > 0))
  expect_true(all(y5e$cum_vaccinated[-1, 2, 3] > 0))
  expect_true(all(y5e$cum_vaccinated[1, 2, 2] == 0))
  expect_true(all(y5e$N[-1, 2, 2] > 0))
})
