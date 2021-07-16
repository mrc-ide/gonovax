devtools::load_all()
tt <- seq(0, 20)
gp <- gono_params(1:500)
y00 <- run_onevax_xvw(c(0, 50), gp, vea = 0, dur = 1e3)
ip <- lapply(y00, restart_params)
y0 <- run_onevax_xvw(tt, gp, ip, vea = 0, dur = 1e3)
y <- run_onevax_xvw(tt, gp, ip, vea = 0.5, uptake = 0.33, strategy = "VoD", dur = 10)

n <- length(tt)

for (i in seq_along(y)) {
  # calculate the baseline diag_a each year
  y0[[i]]$diag_a <- apply(y0[[i]]$cum_diag_a, c(2, 3), diff)
  # calculate the baseline diag_a per person each year
  mean_pop <- (y0[[i]]$N[-1, , "X"] + y0[[i]]$N[-n, , "X"]) / 2
  eta <- rep(c(gp[[i]]$eta_l_t[1], gp[[i]]$eta_h_t[1]), each = n)
  diag_a_per_capita <- y0[[i]]$A[, , "X"] * eta / y0[[i]]$N[, , "X"]
  # multiply the per-capita rate of diag_a to the population structure of the
  # vaccine run to work out the equivalent number of diagnoses in each
  # vaccine strata
  y0[[i]]$diag_a_weighted <- c(diag_a_per_capita) * y[[i]]$N[ , , ]
  y[[i]]$diag_a_weighted <- y[[i]]$A * eta
}

mean_ci <- function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975)))
baseline_diag_a_x <- aggregate(y0, "diag_a_weighted", stratum = "X", f = mean_ci)
vax_diag_a_x <- aggregate(y, "diag_a_weighted", stratum = "X", f = mean_ci)

baseline_diag_a_vw <- aggregate(y0, "diag_a_weighted", stratum = c("V", "W"), f = mean_ci)
vax_diag_a_vw <- aggregate(y, "diag_a_weighted", stratum = c("V", "W"), f = mean_ci)

# baseline_diag_s_x <- aggregate(y0, "diag_s_weighted", stratum = "X")
# vax_diag_s_x <- aggregate(y, "cum_diag_s", as_incid = TRUE, stratum = "X")
# 
# baseline_diag_s_vw <- aggregate(y0, "diag_s_weighted", stratum = c("V", "W"))
# vax_diag_s_vw <- aggregate(y, "cum_diag_s", as_incid = TRUE, stratum = c("V", "W"))

par(mfrow = c(1, 3))
lty <- c(1, 2, 2)
ylim <- c(0, 1e4)
cols <- c("black", "red")
matplot(t(baseline_diag_a_x), type = "l", lty = lty, col = cols[1],
        ylab = "asymp cases", main = 'unvaccinated', ylim = ylim)
matlines(t(vax_diag_a_x), lty = lty, col = cols[2])
legend("topright", fill = cols,
       legend = c("baseline", "vaccine run"))
matplot(t(baseline_diag_a_vw), type = "l", lty = lty, col = cols[1],
        ylab = "asymp cases", main = "ever vaccinated", ylim = ylim)
matlines(t(vax_diag_a_vw), lty = lty, col = cols[2])

cols <- c("blue", "green")
plot(baseline_diag_a_vw["mean", ] - vax_diag_a_vw["mean", ], type = "l", col = cols[1],
     ylab = 'asymp cases averted', ylim = c(0, 7e3), lty = 1)
lines(baseline_diag_a_x["mean", ] - vax_diag_a_x["mean", ], col = cols[2], lty = 1)
legend("topleft", fill = cols,
       legend = c("ever vaccinated", "unvaccinated"))


