test_that("can retrieve data", {
  data <- read_csv(gonovax_file("extdata/observed_data.csv"))
  expect_equal(gonovax_data(), data)
})
