#' Generates a heatmap to illustrate the abundance of organisms across sample microbiomes.
#'
#' `create_heatmap` takes in an abundance matrix to generate a heatmap
#' to illustrate abundance of 16S sequence variants across samples.
#'
#' @param table abundance table generated from `create_abundance_table`
#' @param proportion proportion of table between 0 and 1 to visualize. If 1, will visualize all columns, if 0.5 will visualize the first half of the columns of the table.
#'
#' @returns A ggplot2 heatmap object
#'
#' @examples
#' plots_dir = "data/graphics"
#' \dontrun{
#' heatmap <- create_heatmap(table, plots_dir, proportion = 0.025)
#' }
#'
#' @author {Koji Wong, \email{koji.wong@mail.utoronto.ca}}
#'
#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#'
#' Wickham H, Vaughan D, Girlich M (2023). tidyr: Tidy Messy Data.
#' R package version 1.3.0, <https://CRAN.R-project.org/package=tidyr>.
#'
#' Müller K, Wickham H (2023). tibble: Simple Data Frames.
#' R package version 3.2.1, <https://CRAN.R-project.org/package=tibble>.
#'
#' @import ggplot2
#' @import tibble
#' @import tidyr
#' @export
create_heatmap <- function(table, proportion = 1, low_col = "#F2E7C9", high_col ="#801A86") {
  # Check table class type
  if (!is(table, "matrix")) {
    message = "table parameter must be a 2D matrix where each row is a sample and each column is a 16S sequence."
    stop(message)
  }

  if (proportion <= 0) {
    message = paste("Proportion value must be greater than 0 but received value", proportion)
    stop(message)
  }
  if (proportion > 1) {
    message = paste("Proportion value must be less than or equal to 1 but received value", proportion)
    stop(message)
  }

  len = round(proportion * length(table[1, ]))
  prop_table = table[, 1:len]
  colnames(prop_table) <- paste0("S", 1:len, sep = "")  # Changed column names

  # Clean up our table using tibble and tidyr
  tidy_data <- as.data.frame(prop_table) %>%
    tibble::rownames_to_column(var = "Sample") %>%
    tidyr::gather(key = "Species", value = "Abundance", -Sample)
  # Generate our heatmap using ggplot2
  heatmap <- ggplot2::ggplot(tidy_data, aes(x = Species, y = Sample, fill = Abundance)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = low_col, high = high_col) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Abundance Heatmap")
  print(heatmap)
  return(heatmap)
}
# [END]
