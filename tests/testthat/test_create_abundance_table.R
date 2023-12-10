library("CompHeatmaps")
library("testthat")

testthat::test_that("Can create Abundance Table successfully", {
  result <- readRDS("data/outputs/result")
  table <- CompHeatmaps::create_abundance_table(result)
  # expect that we have a 2D matrix
  expect_equal(length(dim(table)), 2)
  # expect that the number of sequence variants after filtering is 432
  expect_equal(length(table[1, ]) == 432)
  # expect that the number of input samples is 4
  expect_equal(length(table[, 1]) == 4)
})
