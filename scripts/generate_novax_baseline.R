# run no-vaccination model for 20 years from equilibrium
res <- gonovax::run_novax(tt = seq(0, 30), equilib = TRUE)

# extract cumulative incidence (aggregated over activity group and vax status)
# this returns a time x parameter set matrix

cum_incid  <- t(aggregate(res, "cum_incid"))
cum_diag_a <- t(aggregate(res, "cum_diag_a"))
cum_diag_s <- t(aggregate(res, "cum_diag_s"))
screened <- t(aggregate(res, "cum_screened"))

# extract incidence from cumulative incidence and remove cumulative incidence at
# time zero
novax_baseline <- list(incid = apply(cum_incid, 2, diff),
                       cum_incid = cum_incid[-1, ],
                       cum_diag_a = cum_diag_a[-1, ],
                       cum_diag_s = cum_diag_s[-1, ],
                       cum_screened = cum_screened[-1, ])

saveRDS(novax_baseline, "inst/extdata/novax_baseline.rds")
