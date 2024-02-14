##' @name stratum_index_xvwr
##' @title Generate the indices of all xvw trial strata
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
##' @return A list of strata with their indices
##' @export

stratum_index_xvwr <- function(n_diag_rec = 1) {

  # for an n_erlang of 3, and n_diag_rec of 2, the list of indexes returned
  # will be in the following order, where roman numerals refer to n_diag_rec,
  # and arabic numerals refer to erlang:
  # X.I, X.II, V1.I, V1.II, V2.I, V2.II, V3.I, V3.II, W.I, W.II

  ret <- list(X = seq_len(n_diag_rec))

  ret$V <- max(ret$X) + seq_len(n_diag_rec)
  ret$W <- max(ret$V) + seq_len(n_diag_rec)
  ret$R <- max(ret$W) + seq_len(n_diag_rec)

  ret$n_vax <- max(ret$R)
  ret$n_diag_rec <- n_diag_rec

  ret
}

##' @name stratum_index_xvwrh
##' @title Generate the indices of all xvw trial strata
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer giving the number of each X, V(erlang), and W
##' stratum, allowing tracking of diagnosis history. e.g for a n_diag_rec = 2
##' and erlang = 1, there will be X.I, X.II, V1.I, V1.II, W.I, W.II strata.
##' Where '.I' corresponds to never-diagnosed individuals and '.II' is for
##' individuals diagnosed at least once.
##' @return A list of strata with their indices
##' @export
stratum_index_xvwrh <- function(n_diag_rec = 1) {

  # for an n_erlang of 3, and n_diag_rec of 2, the list of indexes returned
  # will be in the following order, where roman numerals refer to n_diag_rec,
  # and arabic numerals refer to erlang:
  # X.I, X.II, V1.I, V1.II, V2.I, V2.II, V3.I, V3.II, W.I, W.II

  ret <- list(X = seq_len(n_diag_rec))

  ret$V <- max(ret$X) + seq_len(n_diag_rec)
  ret$W <- max(ret$V) + seq_len(n_diag_rec)
  ret$R <- max(ret$W) + seq_len(n_diag_rec)
  ret$H <- max(ret$R) + seq_len(n_diag_rec)

  ret$n_vax <- max(ret$H)
  ret$n_diag_rec <- n_diag_rec

  ret
}
