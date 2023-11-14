if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::version()
BiocManager::install("dada2")
BiocManager::install("ShortRead")
BiocManager::install("rhdf5filters")
BiocManager::install("phyloseq")
install.packages("pheatmap")
install.packages("here")

library("here")
library("dada2")
library("phyloseq")
library("pheatmap")
library("ShortRead")

current_dir = getwd()
path <- file.path(current_dir, "data/MiSeq_SOP") # concatenate wd and data files
path
list.files(path) # make sure the fastq files are listed
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
fnFs
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
fout_path <- file.path(current_dir, "/data/output/")

preprocess_16s_data <- function(input_fastq_files, output_files, verbose = FALSE, multithread = FALSE) {
  output_path <- file.path(current_dir, "data/filtered_reads")
  # Quality filter and denoise
  if (verbose) {
    print("Filtering Paired Reads in 16S Data, this may take a couple of minutes...")
  }
  out <- dada2::filterAndTrim(
    sample_list,
    sample_output,
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
    print("Filtering Complete.")
    print("Learning the Error Rates, this may take a couple of minutes...")
  }
  # Learn the Error Rates
  err <- dada2::learnErrors(sample_output_path, multithread = multithread)
  if (verbose) {
    print("Error Learning Complete.")
    print("Piping through dada, this may take a couple of minutes...")
  }
  dada_out <- dada2::dada(sample_output_path, err = err, multithread = multithread)
  return(out)
}


create_abundance_table <- function(filt) {
  # Create a sequence table
  seq_table <- makeSequenceTable(filt)

  # Normalize to relative abundances
  rel_abundance_table <- seq_table / rowSums(seq_table)

  # Optionally, round the values for clarity
  rel_abundance_table <- round(rel_abundance_table, 6)

  return(rel_abundance_table)
}


combine_data <- function(tables) {
  # combine vectors to make a combined abundance matrix
  combined_matrix = cbind(tables)
}

create_heatmap <- function(table) {

}


#===== Examples for Vignettes ========
sample_data_path <- file.path(current_dir, "data/sample_raw_16S_data")
sample_list <-sort(list.files(sample_data_path, pattern=".fastq.gz", full.names = TRUE))
sample_data.names <- sapply(strsplit(basename(sample_list), "_"), `[`, 1)
sample_output_path <- file.path(current_dir, "data/filtered_reads")
sample_output <- file.path(sample_output_path, paste0(sample_data.names, "_filt.fastq.gz"))
names(sample_output) <- sample_data.names

single_entry <- file.path(current_dir, "data/sample_raw_16S_data/HSM5MD4R_P.fastq")
single_entry_output <- file.path(current_dir, "data/filtered_reads/HSM5MD4R_filt.fastq.gz")
table <- preprocess_16s_data(sample_list, sample_output, verbose = TRUE, multithread = TRUE)
out
print(table)
makeSequenceTable(output_files)
plotErrors(err, nominalQ=TRUE)
dada_out <- dada(sample_output_path, err = err)
dada_out

# ===== DADA2 pipeline =====

#Place filtered files in filtered/ subdirectory
sample.names
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
filtFs
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
# Learn the error rates produced from filtering
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Input to dada2 pipeline
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Make sequence table
seqtab <- makeSequenceTable(mergers)
seqtab
dim(seqtab)
# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


# ===== Package functions ======

create_output_dir <- function()