##' @name create_uptake_map_xvw
##' @title Creates uptake mapping array with dimensions n_group x n_vax x n_vax
##' and assigns the relevant primary uptake and booster uptake values defined
##' by the user.
##' @param n_group scalar indicating number of activity groups
##' @param n_vax scalar indicating number the number of stratum in the model
##' @param uptake proportion of the unvaccinated population who accept
##' vaccination
##' @return an array of the uptakes with dimensions n_group x n_vax x n_vax

create_uptake_map_xvw <- function(n_group, n_vax, uptake, idx,
                                  n_diag_rec = 1, screening_or_diagnosis) {
  
  # set up uptake matrix rows = groups, columns = vaccine strata
  u <- array(0, dim = c(n_group, n_vax, n_vax))
  
  for (i in 1: n_diag_rec){
    
    if (screening_or_diagnosis == "screening") {
      temp <- i
    } else if (screening_or_diagnosis == "diagnosis") {
      
      if (i < n_diag_rec) {
        temp <- i + 1
      } else {
        temp <- i
      }
    } else {
      print("uptake map type not specified.")
    }
    
    u[, i, i] <- uptake
    u[, idx$V[temp], i] <- uptake
    
  }
  
  u
}

##' @name create_uptake_map_xvwv
##' @title Creates uptake mapping array with dimensions n_group x n_vax x n_vax
##' and assigns the relevant primary uptake and booster uptake values defined
##' by the user.
##' @param n_group scalar indicating number of activity groups
##' @param n_vax scalar indicating number the number of stratum in the model
##' @param uptake proportion of the unvaccinated population who accept
##' vaccination
##' @return an array of the uptakes with dimensions n_group x n_vax x n_vax

create_uptake_map_xvwv <- function(n_group, n_vax, primary_uptake, booster_uptake,
                                  idx, n_diag_rec = 1, screening_or_diagnosis) {
  
  # set up uptake matrix rows = groups, columns = vaccine strata
  u <- array(0, dim = c(n_group, n_vax, n_vax))
  
  for (i in 1:n_diag_rec){
    
    if (screening_or_diagnosis == "screening") {
      temp <- i
    } else if (screening_or_diagnosis == "diagnosis") {
      
      if (i < n_diag_rec) {
        temp <- i + 1
      } else {
        temp <- i
      }
    } else {
      print("uptake map type not specified.")
    }
    
    u[, i, i] <- primary_uptake
    u[, idx$V[temp], i] <- primary_uptake
    
    u[, idx$W[i], idx$W[i]] <- booster_uptake
    u[, idx$V[temp], idx$W[i]] <- booster_uptake
    
    
  }
  
  u
}

##' @name create_uptake_map_xvwr
##' @title Creates uptake mapping array with dimensions n_group x n_vax x n_vax
##' and assigns the relevant primary uptake and booster uptake values defined
##' by the user.
##' @param n_group scalar indicating number of activity groups
##' @param n_vax scalar indicating number the number of stratum in the model
##' @param uptake proportion of the unvaccinated population who accept
##' vaccination
##' @return an array of the uptakes with dimensions n_group x n_vax x n_vax

create_uptake_map_xvwr <- function(n_group, n_vax, primary_uptake, booster_uptake,
                                   idx, n_diag_rec = 1, screening_or_diagnosis) {
  
  # set up uptake matrix rows = groups, columns = vaccine strata
  u <- array(0, dim = c(n_group, n_vax, n_vax))
  
  for (i in 1: n_diag_rec){
    
    if (screening_or_diagnosis == "screening") {
      temp <- i
    } else if (screening_or_diagnosis == "diagnosis") {
      
      if (i < n_diag_rec) {
        temp <- i + 1
      } else {
        temp <- i
      }
    } else {
      print("uptake map type not specified.")
    }
    
    u[, i, i] <- primary_uptake
    u[, idx$V[temp], i] <- primary_uptake
    
    u[, idx$W[i], idx$W[i]] <- booster_uptake
    u[, idx$R[temp], idx$W[i]] <- booster_uptake
    
    
  }
  
  u
}



##' @name create_uptake_map_xpvwrh
##' @title Creates uptake mapping for the branching XPVWRH model where
##' individuals can move from unvaccinated (X) to vaccinated (V) or partially
##' vaccinated (P) as well as revaccinated from waned (W) to (R) and, and
##' partially vaccinated (P) to fully vaccianted (V). The former
##' reflects the specific indices which are chosen for assigning uptakes.
##' @param array a vaccine map array of dimensions n_group by n_vax by n_vax
##' generated through create_vax_map_branching()
##' @param r1 proportion of population offered vaccine only accepting the first
##' dose, becoming partially vaccinated
##' @param r2 proportion of the population who accepted the first dose of the
##' vaccine who go on to accept the second dose, becoming fully vaccinated
##' @param booster_uptake proportion of the formerly fully vaccinated, waned
##' population who accept a booster vaccination dose
##' @param r2_p proportion of partially vaccinated individuals who receive
##' a second dose when returning to the clinic due to screening or illness
##' @param idx list containing indices of all X, P, V, W, R & H strata and n_vax
##' through vaccine-protected strata until that protection has waned
##' @param n_diag_rec
##' @return an array of the uptakes of same dimensions

create_uptake_map_xpvwrh <- function(array, r1, r2, r2_p, booster_uptake,
                                     idx, n_diag_rec = 1, 
                                     screening_or_diagnosis) {
  
  for (i in 1:n_diag_rec) {
    
    # note, these indices are specific to the branching pattern of xpvwrh
    ## individuals in X accept vaccination of the 1st dose at an uptake of r1
    
    if (screening_or_diagnosis == "screening") {
      temp <- i
    } else if (screening_or_diagnosis == "diagnosis") {
      
      if (i < n_diag_rec) {
        temp <- i + 1
      } else {
        temp <- i
      }
    } else {
      print("uptake map type not specified.")
    }
    
    array[, i, i] <- array[, i, i] * r1
    
    ## individuals entering V (fully vaccinated) also then accept
    ## vaccination with the 2nd dose at an uptake of r2 so the proportion fully
    ## vaccinated is given by r1 * r2
    ## idx$V[1] gives index of the top of the V erlang stack
    array[, idx$V[temp], i] <- array[, idx$V[temp], i] * (r1 * r2)
    
    ## individuals entering P (partially vaccinated) do not then accept
    ## vaccination with the 2nd dose so proportion partially vaccinated is
    ## given by r1 * (1 - r2), where 1 - r2 is the proportion not accepting the
    ## 2nd dose given they have recieved the 1st dose
    ## idx$P[1] gives index of the top of the P erlang stack
    array[, idx$P[temp], i] <- array[, idx$P[temp], i] * (r1 * (1 - r2))
  }
  
  
  ## individuals with only the 1st dose can later accept vaccination with the
  ## 2nd dose at an uptake of r2_p
  ## idx$P gives indices for all P erlang strata, r2_p applies to all equally
  array[, , idx$P] <- array[, , idx$P] * r2_p
  
  
  ## individuals who were fully vaccinated and whose immunity has waned (W)
  ## can accept vaccination with a single booster dose at an uptake of
  ## booster_uptake
  ## idx$W gives the the index for (W)
  
  array[, , idx$W] <- array[, , idx$W] * booster_uptake
  
  
  # values must be positive - otherwise negative values in this array will
  # cancel those in the vos and vod arrays = incorrect vaccination
  abs(array)
}

