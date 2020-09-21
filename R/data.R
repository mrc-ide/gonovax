##' Get annual time series of GUMCAD diagnoses and testing figures
##'
##' @title Get GUMCAD data
##' @return A data.frame containing the gonovax year (i.e. year after 2007),
##' the year, the total diagnoses in MSM, the total tests in MSM and the source
##' of the data
##'
##' @export
gumcad_data <- function() {
if (is.null(cache$data)) {
  cache$data <- read_csv(gonovax_file("extdata/observed_gumcad.csv"))
}
  cache$data
}
