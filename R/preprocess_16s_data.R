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
#' \dontrun{
#' library("CompHeatmaps")
#' input <- system.file("extdata/sample_raw_16S_data", package = "CompHeatmaps")
#' output <- system.file("extdata/filtered_reads", package = "CompHeatmaps")
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
  return(dada_out)
}
# [END]
