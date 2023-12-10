library("CompHeatmaps")
library("testthat")

testthat::test_that("Preprocess works successfully", {

  input_dir = "../../inst/extdata/sample_raw_16S_data"
  output_dir = "../../inst/extdata/filtered_reads"
  result <- CompHeatmaps::preprocess_16s_data(input_dir,
                                              output_dir,
                                              verbose=TRUE,
                                              multithread=TRUE)

  expect_error(result<- CompHeatmaps::preprocess_16s_data("/invalid/dir", output_dir))
  expect_error(result<- CompHeatmaps::preprocess_16s_data(input_dir, "/invalid/dir"))
})
