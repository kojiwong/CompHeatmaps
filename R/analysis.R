if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("phyloseq")
library("dada2")
packageVersion("dada2")
path <- paste(getwd(), "/data/MiSeq_SOP", sep="") # concatenate wd and data files
list.files(path) # make sure the fastq files are listed
