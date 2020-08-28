tt <- c(0, 30)
novax <- gonovax::run_novax(tt = tt)
saveRDS(novax, "inst/extdata/novax_equilib.rds")
