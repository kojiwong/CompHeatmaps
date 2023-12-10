#' Preprocess 16S Sequence Data using Dada2 Package
#'
#' `preprocess_16s_data` preprocesses 16S rRNA single or paired-end reads
#' that are in FASTQ format. This function calls the `filterAndTrim` function from
#' the dada2 package which filters and trims fastq file inputs.
#'
#' @param input_dir directory containing desired input files. Must be in fastq.gz, fastq.bz2, or .fastq file
#' @param output_dir target directory where preprocessed data will be saved. This is the same directory that will be used to learn error rates.
#' @param verbose set to true for explicit updates during preprocessing
#' @param multithread set to true if you wish to utilize multiple CPU cores for extra speed.
#'
#' @returns a dada output file
#'
#' @examples
#' # Example:
#' # Directory containing raw single-end 16S reads
#' input <- system.file("inst/extdata/sample_raw_16S_data", package = "CompHeatmaps")
#' output <- system.file("/inst/extdata/filtered_reads", package = "CompHeatmaps")
#' \dontrun{
#' library("CompHeatmaps")
#' # Set verbose to TRUE for updates on progress in terminal and multithread to TRUE for faster computation if you have multiple cores.
#' CompHeatmaps::preprocess_16s_data(input_dir = input, output_dir = output, verbose = TRUE, multithread = TRUE)
#' length(result) # 4
#' class(result) # list
#' class(result[[1]]) # dada2
#' }
#'
#' @author
#' Koji Wong
#'
#' @import dada2
#' @export
preprocess_16s_data <- function(input_dir, output_dir, verbose = FALSE, multithread = FALSE) {
  # Check if directories are valid
  if (!is.character(input_dir) || length(input_dir) == 0 || !dir.exists(input_dir)) {
    stop("Specified directory containing input data is invalid. Please provide a valid and existing directory path.")
  }
  if (!is.character(output_dir) || length(output_dir) == 0 || !dir.exists(output_dir)) {
    stop("Specified output directory is invalid. Please provide a valid and existing directory path.")
  }
  if (input_dir == output_dir) {
    stop("Output directory must be distinct from input directory.")
  }
  # Check if input directory contains valid files
  files <- list.files(input_dir)
  valid_formats <- c(".fastq.gz", ".fastq.bz2", ".fastq")
  # Check if there are files with the specified formats
  matching_files <- grep(paste(valid_formats, collapse = "|"), files, value = TRUE)
  if (length(matching_files) == 0) {
    stop("No files with valid formats (.fastq.gz, .fastq.bz2, .fastq) found in the directory.")
  }

  # Quality filter and denoise
  if (verbose) {
    print("Filtering Paired Reads in 16S Data, this may take a couple of minutes...")
  }
  # Key step. Calls the dada2 filterAndTrim function which takes the raw 16S data
  # from our input directory and filters and trims each fastq file. Outputs compressed
  # fastq files containing trimmed reads which passed the filters.
  out <- dada2::filterAndTrim(
    input_dir,
    output_dir,
    trimLeft = 0,
    truncLen = 140,
    maxN = 0,
    maxEE = 2,
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = multithread
  )
  if (verbose) {
    print("Filtering Step complete.")
    print("Learning the Error Rates, this may take a couple of minutes...")
  }
  # Learn the Error Rates
  err <- dada2::learnErrors(output_dir, multithread = multithread)
  if (verbose) {
    print("Error Learning Step complete.")
    print("Piping through dada, this may take a couple of minutes...")
  }
  # Save our dada2 output and return it
  dada_out <- dada2::dada(output_dir, err = err, multithread = multithread)
  if (verbose) {
    print("Preprocessing of 16S data done.")
  }
  saveRDS(dada_out, system.file("data/precomputed/result", package = "CompHeatmaps"))
  return(dada_out)
}

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
#' /dontrun{
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

