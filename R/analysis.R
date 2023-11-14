if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::version()
BiocManager::install("phyloseq")
BiocManager::install("dada2")
BiocManager::install("ShortRead")
# BiocManager::install("ComplexHeatmap")
# install.packages("pheatmap")
install.packages("here")
install.packages(c("circlize", "ComplexUpset"))
install.packages("gplots")

# Package developing tools
library("devtools")
library("usethis")

# Packages used by library
library("circlize")
library("ComplexUpset")
library("here")
library("dada2")
library("phyloseq")
# library("pheatmap")
library("ShortRead")
library("gplots")

# Set current_dir to refer throughout package
current_dir = getwd()
print(current_dir)

# ===== Package functions ======
#' Preprocess 16S Data
#'
#' `preprocess_16s_data` preprocesses 16S rRNA single or paired-end reads
#' that are in FASTQ format.
#'
#' @param input_fastq_files directory containing desired input files
#' @param output_files directory where output files will go
#' @param verbose set to true for explicit updates in processing
#' @param multithread set to true if you wish to utilize multiple CPU cores
#' @returns a dada output file
#' @importFrom dada2
#' @export
preprocess_16s_data <- function(inputs, outputs, verbose = FALSE, multithread = FALSE) {
  # Set current_dir to refer throughout package
  current_dir = getwd()
  print(current_dir)

  output_path <- file.path(current_dir, "data/filtered_reads")
  # Quality filter and denoise
  if (verbose) {
    print("Filtering Paired Reads in 16S Data, this may take a couple of minutes...")
  }

  out <- dada2::filterAndTrim(
    inputs,
    outputs,
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
  err <- dada2::learnErrors(sample_output_path, multithread = multithread)
  if (verbose) {
    print("Error Learning Step complete.")
    print("Piping through dada, this may take a couple of minutes...")
  }
  dada_out <- dada2::dada(sample_output_path, err = err, multithread = multithread)
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


# #===== Examples for Vignettes ========
# current_dir <- getwd()
# sample_data_path <- file.path(current_dir, "data/sample_raw_16S_data")
# sample_output_path <- file.path(current_dir, "data/filtered_reads")

# sample_list <-sort(list.files(sample_data_path, pattern=".fastq.gz", full.names = TRUE))
# sample_data.names <- sapply(strsplit(basename(sample_list), "_"), `[`, 1)
# sample_output <- file.path(sample_output_path, paste0(sample_data.names, "_filt.fastq.gz"))
# sample_output
# names(sample_output) <- sample_data.names
# dada_output <- preprocess_16s_data(sample_list, sample_output, verbose = TRUE, multithread = TRUE)
# dada_output

# # Make a sequence table
# seq_table <- create_abundance_table(dada_output)
# rownames(seq_table) <- sample_data.names
# seq_table
# dim(seq_table)
# table(nchar(getSequences(seq_table)))

# # Remove chimeras
# seq_table.nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=TRUE, verbose=TRUE)
# dim(seq_table.nochim)

# # Assign taxonomy
# taxon_ref_path = "data/ref_datasets/silva_nr99_v138.1_train_set.fa.gz"
# taxon_database_path <- file.path(current_dir, taxon_ref_path)
# taxon_database_path
# taxa <- assignTaxonomy(seq_table.nochim, taxon_database_path)
# taxa.print <- taxa
# row_heights <- rep(0.2, nrow(seq_table))
# rownames(taxa.print) <- NULL
# head(taxa.print)
# taxa.print

# # Plotting


# ===== Developing Package =====
usethis::use_mit_license("Koji Wong")
usethis::use_package("dada2", type = "Imports", min_version = "1.30")
usethis::use_package("gplots", type = "Imports", min_version = "3.1.3")
usethis::use_package("here", type = "Imports", min_version = "1.0.1")
usethis::use_package("circlize", type = "Imports", min_version = "0.4.15")
usethis::use_package("phyloseq", type = "Imports", min_version = "1.46.0")
usethis::use_package("ComplexUpset", type = "Imports", min_version = "1.3.3")
usethis::use_package("ShortRead", type = "Imports", min_version = "1.60.0")