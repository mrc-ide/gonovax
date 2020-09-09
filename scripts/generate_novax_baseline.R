# run no-vaccination model for 20 years from equilibrium
res <- gonovax::run_novax(tt = seq(0, 30), equilib = TRUE)

# extract cumulative and incident flows
novax_baseline <- extract_flows(res)

# save results
saveRDS(novax_baseline, "inst/extdata/novax_baseline.rds")
