

vax_params1 <- function(dur = 4, ve = 1, vs = 0, vd = 0, eff = 0.31) {
  z <- 1 / dur
  list(n_vax = 2,
       ve    = c(1 - ve, ve),
       vs    = vs,
       vd    = vd,
       eff   = c(0, eff),
       v     = rbind(c(1, 0),
                     c(-1, 0)),
       w     = rbind(c(0,  z),
                     c(0, -z))
  )
}