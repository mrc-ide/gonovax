stratum_index_xvw <- function(n_erlang = 1, n_diag_rec = 1, strategy = NULL) {
  
  # for an n_erlang of 3, and n_diag_rec of 2, the list of indexes returned
  # will be in the following order, where roman numerals refer to n_diag_rec,
  # and arabic numerals refer to erlang:
  # X.I, X.II, V1.I, V1.II, V2.I, V2.II, V3.I, V3.II, W.I, W.II
  
  ret <- list(X = seq_len(n_diag_rec))
  ret$V <- max(ret$X) + seq_len(n_erlang * n_diag_rec)
  ret$W <- max(ret$V) + seq_len(n_diag_rec)
  ret$n_vax <- max(ret$W)
  n_vax <- ret$n_vax
  
  # strata people are diagnosed from
  ret$diagnosedfrom = seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 0]
  
  # strata people are diagnosed to
  ret$diagnosedto = seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 1]
  
  # strata people are vaccinated (before entry) from and to
  ret$vaccinatedfrom_vbe = ret$X
  ret$vaccinatedto_vbe = ret$V
  
  # strata people are vaccinated (on screening) from and to
  if (!is.null(strategy)){
  if (!is.null(strategy) & (strategy == "VaH" | strategy ==  "VaHonly")){
    ret$vaccinatedfrom_vos = ret$X[-1]
    ret$vaccinatedto_vos = ret$V[2:n_diag_rec]
  } else {
    ret$vaccinatedfrom_vos = ret$X
    ret$vaccinatedto_vos = ret$V
  }
  }
  
  #strata people are vaccinated (on diagnosis) from and to
  ret$vaccinatedfrom_vod = ret$X
  ret$vaccinatedto_vod = ret$V
  ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] = 
    ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] + 1
    
  ret
}

stratum_index_xvwv <- function(n_erlang = 1, n_diag_rec = 1, strategy = NULL) {
  
  # for an n_erlang of 3, and n_diag_rec of 2, the list of indexes returned
  # will be in the following order, where roman numerals refer to n_diag_rec,
  # and arabic numerals refer to erlang:
  # X.I, X.II, V1.I, V1.II, V2.I, V2.II, V3.I, V3.II, W.I, W.II
  
  ret <- list(X = seq_len(n_diag_rec))
  ret$V <- max(ret$X) + seq_len(n_erlang * n_diag_rec)
  ret$W <- max(ret$V) + seq_len(n_diag_rec)
  ret$n_vax <- max(ret$W)
  n_vax <- ret$n_vax
  
  # strata people are diagnosed from
  ret$diagnosedfrom = seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 0]
  
  # strata people are diagnosed to
  ret$diagnosedto = seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 1]
  
  # strata people are vaccinated (before entry) from and to
  ret$vaccinatedfrom_vbe = ret$X
  ret$vaccinatedto_vbe = ret$V
  
  # strata people are vaccinated (on screening) from and to
  if (!is.null(strategy)){
    if (!is.null(strategy) & (strategy == "VaH" | strategy ==  "VaHonly")){
      ret$vaccinatedfrom_vos = c(ret$X[-1], ret$W[-1])
      ret$vaccinatedto_vos = c(ret$V[2:n_diag_rec], ret$V[2:n_diag_rec])
    } else {
      ret$vaccinatedfrom_vos = c(ret$X, ret$W)
      ret$vaccinatedto_vos = c(ret$V, ret$V)
    }
  }
  
  #strata people are vaccinated (on diagnosis) from and to
  ret$vaccinatedfrom_vod = c(ret$X, ret$W)
  ret$vaccinatedto_vod = c(ret$V, ret$V)
  ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] = 
    ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] + 1
  
  ret
}







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
