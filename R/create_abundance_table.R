#' Creates a matrix table of abundance counts of 16S reads across samples.
#'
#' `create_abundance_table` is takes in as input a list of dada2 objects
#' that have been filtered to create an abundance matrix table with dada2
#' which will be used for heatmap generation. Performs some normalization
#' and filtering to simplify results. Removes sequence variants identified as
#' bimeric from the table to clean up further.
#' @param filtered_results the filtered results produced by dada2, a list of dada objects
#' @param multithread set to TRUE to compute with multithreading for faster generation, only works with multiple cores
#' @returns an abundance table
#'
#' @examples
#' library("CompHeatmaps")
#' \dontrun{
#' # get filtered_reads from running `preprocess_16s_data()` on dataset
#' table <- CompHeatmaps::create_abundance_table(filtered_reads)
#' table[[1]] # -0.632696, abundance count after normalized and scaled
#' dim(table) # 4 432, 4 rows each a sample in our input directory and 432 columns each a 16S sequence
#' # can further explore data with `table[1, ]`, `table[, 1]`
#' }
#'
#' @author {Koji Wong}
#'
#' @import dada2
#' @export
create_abundance_table <- function(filtered_results, multithread = FALSE) {
  # type checking for creating abundance table
  if (class(filtered_results) != "list") {
    stop("filtered_results must be of type list. Please use the output from preprocess_16s_data() to produce input for this function.")
  }

  for (i in 1:length(filtered_results)) {
    if (class(filtered_results[[i]]) != "dada") {
      message = paste("Wrong parameter type. filtered_results must be a list consisting of dada objects to create an abundance table from but element", i, "of filtered_results has type", class(filtered_results[[i]]))
      stop(message)
    }
  }

  # Create a sequence table
  seq_table <- dada2::makeSequenceTable(filtered_results, orderBy = "abundance")

  # Keep only abundances above threshold of 10
  seq_table <- seq_table[, rowSums(seq_table) > 10]

  # Remove chimeras
  seq_table_nochim <- dada2::removeBimeraDenovo(seq_table, method = "consensus", multithread = FALSE, verbose = FALSE)

  # Normalize to relative abundances
  rel_abundance_table <- scale(seq_table_nochim)

  # Optionally, round the values for clarity
  rel_abundance_table <- round(rel_abundance_table, 6)

  return(rel_abundance_table)
}
# [END]
