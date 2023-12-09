# ===== Package functions ======
#' Preprocess 16S Data
#'
#' `preprocess_16s_data` preprocesses 16S rRNA single or paired-end reads
#' that are in FASTQ format.
#'
#' @param input_dir directory containing desired input files. Must be in fastq.gz, fastq.bz2, or .fastq file
#' @param output_dir target directory where preprocessed data will be saved. This is the same directory that will be used to learn error rates.
#' @param verbose set to true for explicit updates during preprocessing
#' @param multithread set to true if you wish to utilize multiple CPU cores
#' @returns a dada output file
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
  return(dada_out)
}

#' Creates an abundance matrix table
#'
#' `create_abundance_table` is used to create an abundance matrix table with dada2
#' which will be used for heatmap generation. Performs some normalization
#' and filtering to simplify results.
#'
#' @param filtered_results the filtered results produced by dada2, a list of dada objects
#' @returns an abundance table
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
#' `create_heatmap` takes in an abundance table to generate a heatmap
#' to illustrate abundance across samples.
#' @param table abundance table generated from `create_abundance_table`
#' @param output_dir directory to save output heatmap object to
#' @param proportion proportion of table between 0 and 1 to visualize. If 1, will visualize all columns, if 0.5 will visualize the first half of the columns of the table.
#' @returns A heatmap object
#' @import gplots
#' @imports ggplot2
#' @export
create_heatmap <- function(table, output_dir, proportion = 1) {
  # Check table class type
  if (!is(table, "matrix")) {
    message = "table parameter must be a 2D matrix where each row is a sample and each column is a 16S sequence."
    stop(message)
  }
  # Check if output directory exists
  if (!is.character(output_dir) || length(output_dir) == 0 || !dir.exists(output_dir)) {
    stop("Specified output directory is invalid. Please provide a valid and existing directory path.")
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

  # Generate our heatmap using gplots heatmap function
  heatmap <- gplots::heatmap.2(as.matrix(table[, 1:len]),
              scale = "row",
              trace = "none",
              main = "Abundance Heatmap",
              key = TRUE,
              keysize = 1,
              density.info = "none",
              margins = c(5, 10),
              Colv = TRUE,
              Rowv = TRUE,
              col = heat.colors(112))
  ggplot2::
  saveRDS(heatmap, paste(output_dir, "heatmap", sep = "/"))
}
