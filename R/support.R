gono_params <- function() {
  if (is.null(cache$gono_params)) {
    cache$gono_params <-
      read_csv(sircovid_file("extdata/gono_params.csv"))
  }
  cache$sgono_params
}
