test_that("can retrieve data", {
  data <- read_csv(gonovax_file("extdata/observed_gumcad.csv"))
  expect_equal(gumcad_data(), data)
})
