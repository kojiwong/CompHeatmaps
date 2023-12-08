# ===== Package functions ======
#' Preprocess 16S Data
#'
#' `preprocess_16s_data` preprocesses 16S rRNA single or paired-end reads
#' that are in FASTQ format.
#'
#' @param input_dir directory containing desired input files. Must be in fastq.gz, fastq.bz2, or .fastq file
#' @param output_dir target directory where preprocessed data will be saved
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
#' @param filtered_results the filtered results produced by dada2
#' @returns an abundance table
#' @import dada2
#' @export
create_abundance_table <- function(filtered_results) {
  # Create a sequence table
  seq_table <- dada2::makeSequenceTable(filtered_results, orderBy = "abundance")

  # Keep only abundances above threshold of 10
  seq_table <- seq_table[, rowSums(seq_table) > 10]

  # Normalize to relative abundances
  rel_abundance_table <- seq_table / rowSums(seq_table)

  # Optionally, round the values for clarity
  rel_abundance_table <- round(rel_abundance_table, 6)

  # Remove chimeras
  seq_table.nochim <- dada2::removeBimeraDenovo(seq_table, method="consensus", multithread=TRUE, verbose=TRUE)
  dim(seq_table.nochim)

  return(seq_table.nochim)
}

#' Generates heatmap
#'
#' `create_heatmap` takes in an abundance table to generate a heatmap
#' to illustrate abundance across samples.
#' @param table abundance table generated from `create_abundance_table`
#' @returns None. Generates a heatmap
#' @import gplots
#' @export
create_heatmap <- function(table) {
  gplots::heatmap.2(as.matrix(table),
          scale = "row",
          trace = "none",
          main = "Abundance Heatmap",
          key = TRUE,
          keysize = 1,
          density.info = "none",
          margins = c(5, 10),
          Colv = TRUE,
          Rowv = TRUE,
          col = heat.colors(256))
}
