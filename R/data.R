##' Get annual time series of GUMCAD diagnoses and testing figures, and annual
##' proportion symptomatic from GRASP
##'
##' @title Get gonovax data
##' @return A data.frame containing the gonovax year (i.e. year after 2007),
##' the year, the total diagnoses in MSM, the total tests in MSM and the source
##' of the data, the proportion of symptomatic diagnoses
##'
##' @export
gonovax_data <- function() {
if (is.null(cache$data)) {
  cache$data <- read_csv(gonovax_file("extdata/observed_data.csv"))
}
  cache$data
}
