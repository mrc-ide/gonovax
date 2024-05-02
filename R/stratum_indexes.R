##' @name stratum_index_xvw
##' @title Generate the indices of all xvwv trial strata
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @param strategy string of vaccination strategy being considered
##' @return A list of strata with their indicies
##' @export
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
  ret$V1 <- ret$V[1:n_diag_rec]

  # strata people are diagnosed from
  ret$diagnosedfrom <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 0]
  
  #new way - by time
  #ret$diagnosedfrom <- seq_len(n_vax)
  
  
  # strata people are diagnosed to
  ret$diagnosedto <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 1]

  #new way - by time
  
  # if (n_diag_rec == 1){
  #   ret$diagnosedto <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 1]
  # } else if (n_diag_rec == 2){
  #   ret$diagnosedto <- rep(seq_len(n_vax)[seq_len(n_vax)%%n_diag_rec == 0], each = n_diag_rec)
  # } else{
  #   ret$diagnosedto <- rep(seq_len(n_vax)[seq_len(n_vax)%%n_diag_rec == 2], each = n_diag_rec)
  # }
  
  
  
  # strata people are vaccinated (before entry) from and to
  ret$vaccinatedfrom_vbe <- ret$X
  ret$vaccinatedto_vbe <- ret$V1

  # strata people are vaccinated (on screening) from and to
  if (!is.null(strategy)) {
    if (!is.null(strategy) && (strategy == "VaH" || strategy ==  "VaHonly" || strategy == "VaH+VoN")) {
      ret$vaccinatedfrom_vos <- ret$X[-1]
      ret$vaccinatedto_vos <- ret$V1[-1]
    } else {
      ret$vaccinatedfrom_vos <- ret$X
      ret$vaccinatedto_vos <- ret$V1
    }
  }

  #strata people are vaccinated (on diagnosis) from and to
  ret$vaccinatedfrom_vod <- ret$X
  ret$vaccinatedto_vod <- ret$V
  ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] <-
    ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] + 1
  
  
  
  ret$vaccinatedfrom_vopn <- ret$X
  ret$vaccinatedto_vopn <- ret$V

  ret
}


##' @name stratum_index_xvwv
##' @title Generate the indices of all xvwv trial strata
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @param strategy string of vaccination strategy being considered
##' @return A list of strata with their indicies
##' @export
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

  ret$V1 <- ret$V[1:n_diag_rec]

  # strata people are diagnosed from
  ret$diagnosedfrom <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 0]

  # strata people are diagnosed to
  ret$diagnosedto <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 1]

  # strata people are vaccinated (before entry) from and to
  ret$vaccinatedfrom_vbe <- ret$X
  ret$vaccinatedto_vbe <- ret$V1

  # strata people are vaccinated (on screening) from and to
  if (!is.null(strategy)) {
    if (!is.null(strategy) && (strategy == "VaH" || strategy ==  "VaHonly"  || strategy == "VaH+VoN")) {
      ret$vaccinatedfrom_vos <- c(ret$X[-1], ret$W[-1])
      ret$vaccinatedto_vos <- c(ret$V1[-1], ret$V1[-1])
    } else {
      ret$vaccinatedfrom_vos <- c(ret$X, ret$W)
      ret$vaccinatedto_vos <- c(ret$V1, ret$V1)
    }
  }

  #strata people are vaccinated (on diagnosis) from and to
  ret$vaccinatedfrom_vod <- c(ret$X, ret$W)
  ret$vaccinatedto_vod <- c(ret$V, ret$V)
  ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] <-
    ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] + 1
  
  ret$vaccinatedfrom_vopn <- c(ret$X, ret$W)
  ret$vaccinatedto_vopn <- c(ret$V, ret$V)
  

  ret
}


##' @name stratum_index_xvwr
##' @title Generate the indices of all xvwr trial strata
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @param strategy string of vaccination strategy being considered
##' @return A list of strata with their indicies
##' @export

stratum_index_xvwr <- function(n_erlang = 1, n_diag_rec = 1, strategy = NULL) {
  # for an n_erlang of 3, and n_diag_rec of 2, the list of indexes returned
  # will be in the following order, where roman numerals refer to n_diag_rec,
  # and arabic numerals refer to erlang:
  # X.I, X.II, V1.I, V1.II, V2.I, V2.II, V3.I, V3.II, W.I, W.II,
  # R1.I, R1.II, R2.I, R2.II, R3.I, R3.II,

  ret <- list(X = seq_len(n_diag_rec))
  ret$V <- max(ret$X) + seq_len(n_erlang * n_diag_rec)
  ret$W <- max(ret$V) + seq_len(n_diag_rec)
  ret$R <- max(ret$W) + seq_len(n_erlang * n_diag_rec)
  ret$n_vax <- max(ret$R)
  n_vax <- ret$n_vax

  ret$V1 <- ret$V[1:n_diag_rec]
  ret$R1 <- ret$R[1:n_diag_rec]

  # strata people are diagnosed from
  ret$diagnosedfrom <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 0]

  # strata people are diagnosed to
  ret$diagnosedto <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 1]

  # strata people are vaccinated (before entry) from and to
  ret$vaccinatedfrom_vbe <- ret$X
  ret$vaccinatedto_vbe <- ret$V1

  # strata people are vaccinated (on screening) from and to
  if (!is.null(strategy)) {
    if (!is.null(strategy) && (strategy == "VaH" || strategy ==  "VaHonly"  || strategy == "VaH+VoN")) {
      ret$vaccinatedfrom_vos <- c(ret$X[-1], ret$W[-1])
      ret$vaccinatedto_vos <- c(ret$V1[-1], ret$R1[-1])
    } else {
      ret$vaccinatedfrom_vos <- c(ret$X, ret$W)
      ret$vaccinatedto_vos <- c(ret$V1, ret$R1)
    }
  }

  #strata people are vaccinated (on diagnosis) from and to
  ret$vaccinatedfrom_vod <- c(ret$X, ret$W)
  ret$vaccinatedto_vod <- c(ret$V, ret$R)
  ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] <-
    ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] + 1
  
  ret$vaccinatedfrom_vopn <- c(ret$X, ret$W)
  ret$vaccinatedto_vopn <- c(ret$V, ret$R)

  ret
}






##' @name stratum_index_xvwrh
##' @title Generate the indices of all xvwrh trial strata
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @param strategy string of vaccination strategy being considered
##' @return A list of strata with their indicies
##' @export
stratum_index_xvwrh <- function(n_erlang = 1, n_diag_rec = 1, strategy = NULL) {
  # for an n_erlang of 3, and n_diag_rec of 2, the list of indexes returned
  # will be in the following order, where roman numerals refer to n_diag_rec,
  # and arabic numerals refer to erlang:
  # X.I, X.II, V1.I, V1.II, V2.I, V2.II, V3.I, V3.II, W.I, W.II,
  # R1.I, R1.II, R2.I, R2.II, R3.I, R3.II, H.I, H.II,

  ret <- list(X = seq_len(n_diag_rec))

  ret$V <- max(ret$X) + seq_len(n_erlang * n_diag_rec)
  ret$W <- max(ret$V) + seq_len(n_diag_rec)
  ret$R <- max(ret$W) + seq_len(n_erlang * n_diag_rec)
  ret$H <- max(ret$R) + seq_len(n_diag_rec)

  ret$V1 <- ret$V[1:n_diag_rec]
  ret$R1 <- ret$R[1:n_diag_rec]

  ret$n_vax <- max(ret$H)
  n_vax <- ret$n_vax

  # strata people are diagnosed from
  ret$diagnosedfrom <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 0]

  # strata people are diagnosed to
  ret$diagnosedto <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 1]

  # strata people are vaccinated (before entry) from and to
  ret$vaccinatedfrom_vbe <- ret$X
  ret$vaccinatedto_vbe <- ret$V1

  # strata people are vaccinated (on screening) from and to
  if (!is.null(strategy)) {
    if (!is.null(strategy) && (strategy == "VaH" || strategy ==  "VaHonly"  || strategy == "VaH+VoN")) {
      ret$vaccinatedfrom_vos <- c(ret$X[-1], ret$W[-1])
      ret$vaccinatedto_vos <- c(ret$V1[-1], ret$R1[-1])
    } else {
      ret$vaccinatedfrom_vos <- c(ret$X, ret$W)
      ret$vaccinatedto_vos <- c(ret$V1, ret$R1)
    }
  }

  #strata people are vaccinated (on diagnosis) from and to
  ret$vaccinatedfrom_vod <- c(ret$X, ret$W)
  ret$vaccinatedto_vod <- c(ret$V, ret$R)
  ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] <-
    ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] + 1
  
  ret$vaccinatedfrom_vopn <- c(ret$X, ret$W)
  ret$vaccinatedto_vopn <- c(ret$V, ret$R)

  ret
}


##' @name stratum_index_xpvwrh
##' @title Generate the indices of all xpvwrh strata
##' @param n_erlang integer giving the number of transitions that need to be
##' made through vaccine-protected strata until that protection has waned
##' @param n_diag_rec integer for the number of diagnosis history substrata
##' @param strategy string of vaccination strategy being considered
##' @return A list of strata with their indicies
##' @export

stratum_index_xpvwrh <- function(n_erlang = 1, n_diag_rec = 1,
                                 strategy = NULL) {
  # for an n_erlang of 3, and n_diag_rec of 2, the list of indexes returned
  # will be in the following order, where roman numerals refer to n_diag_rec,
  # and arabic numerals refer to erlang:
  # X.I, X.II, P1.I, P1.II, P2.I, P2.II, P3.I, P3.II
  # V1.I, V1.II, V2.I, V2.II, V3.I, V3.II, W.I, W.II,
  # R1.I, R1.II, R2.I, R2.II, R3.I, R3.II, H.I, H.II,

  ret <- list(X = 1:n_diag_rec)

  ret$P <- max(ret$X) + seq_len(n_erlang * n_diag_rec)
  ret$V <- max(ret$P) + seq_len(n_erlang * n_diag_rec)
  ret$W <- max(ret$V) + (1:n_diag_rec)
  ret$R <- max(ret$W) + seq_len(n_erlang * n_diag_rec)
  ret$H <- max(ret$R) + (1:n_diag_rec)
  ret$n_vax <- max(ret$H)

  n_vax <- ret$n_vax

  ret$P1 <- ret$P[1:n_diag_rec]
  ret$V1 <- ret$V[1:n_diag_rec]
  ret$R1 <- ret$R[1:n_diag_rec]

  # strata people are diagnosed from
  ret$diagnosedfrom <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 0]

  # strata people are diagnosed to
  ret$diagnosedto <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 1]

  # strata people are vaccinated (before entry) from and to
  ret$vaccinatedfrom_vbe <- c(ret$X, ret$X)
  ret$vaccinatedto_vbe <- c(ret$P1, ret$V1)

  # strata people are vaccinated (on screening) from and to
  if (!is.null(strategy)) {
    if (!is.null(strategy) && (strategy == "VaH" || strategy ==  "VaHonly"  || strategy == "VaH+VoN")) {

      ret$vaccinatedfrom_vos <- c(ret$X[-1], ret$X[-1],
                                  ret$P[- (seq_len(n_erlang * n_diag_rec) %%
                                             n_diag_rec == 1)] , ret$W[-1])
      ret$vaccinatedto_vos <- c(ret$P1[-1], ret$V1[-1],
                                rep(ret$V1[-1], n_erlang), ret$R[-1])
    } else {
      ret$vaccinatedfrom_vos <- c(ret$X, ret$X, ret$P, ret$W)
      # note: V[1] is repeated n_erlang + 1 times as as well as un-vaccinated
      # individuals (X) becoming fully vacc (V), partially vaccinated people
      # in P can become fully vaccinated from any of the P erlang compartments
      ret$vaccinatedto_vos <- c(ret$P1, ret$V1, rep(ret$V1, n_erlang), ret$R)
    }
  }

  #strata people are vaccinated (on diagnosis) from and to
  ret$vaccinatedfrom_vod <- c(ret$X, ret$X, ret$P, ret$W)
  ret$vaccinatedto_vod <- c(ret$P1, ret$V1, rep(ret$V1, n_erlang), ret$R)
  ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] <-
    ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] + 1
  
  ret$vaccinatedfrom_vopn <- c(ret$X, ret$X, ret$P, ret$W)
  ret$vaccinatedto_vopn <- c(ret$P1, ret$V1, rep(ret$V1, n_erlang), ret$R)

  ret
}
