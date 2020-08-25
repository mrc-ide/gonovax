res <- gonovax::run_novax(tt = seq(0, 10), equilib = TRUE)

cum_incid <- sapply(res, function(x) apply(x[["cum_incid"]], 1, sum))

novax_baseline <- list(incid = apply(cum_incid, 2, diff),
                       cum_incid = cum_incid[-1, ])

saveRDS(novax_baseline, "inst/extdata/novax_baseline.rds")
