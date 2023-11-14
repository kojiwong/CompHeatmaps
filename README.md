
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CompHeatmaps

<!-- badges: start -->
<!-- badges: end -->

A large problem in conducting metatranscriptomic analysis on microbial
communities is dealing with the complexity of the large amount of data
we are working with. Heatmaps are a useful visualization that illustrate
the abundance across these samples at a glance. However, getting from
raw sequence reads to heatmaps is quite a cumbersome process.

The goal of CompHeatmaps is to simplify this complex process by
generating heatmaps as easy as 1, 2, 3.

This package was created in R version 4.3.1

## Installation

You can install the development version of CompHeatmaps from
[GitHub](https://github.com/kojiwong/CompHeatmaps) with:

``` r
# install.packages("devtools")
devtools::install_github("kojiwong/CompHeatmaps")
```

## Overview

The main components of this R package include the functions
`preprocess_16s_data`, `create_abundance_table`, and `create_heatmap`

    ls("package:CompHeatmaps")
    data(package = "CompHeatmaps")
    browseVignettes("CompHeatmaps")

## Contributions

## Example

This is a basic example which shows you how to solve a common problem:

# `` {r example} # library(CompHeatmaps) # ## basic example code #  # # In three simple steps, we go from raw 16S rRNA reads to a heatmap illustrating the abundance  # # of short sequence reads across samples. # current_dir = getwd() # # get our paths # sample_data_path <- file.path(current_dir, "data/sample_raw_16S_data") # sample_output_path <- file.path(current_dir, "data/filtered_reads") #  # # get list of files in each path by their extension # sample_list <-sort(list.files(sample_data_path, pattern=".fastq.gz", full.names = TRUE)) # sample_data.names <- sapply(strsplit(basename(sample_list), "_"), `[`, 1) # sample_output <- file.path(sample_output_path, paste0(sample_data.names, "_filt.fastq.gz")) # sample_output # names(sample_output) <- sample_data.names #  # # run through our functions # processed_data <- CompHeatmaps::preprocess_16s_data(sample_list, sample_output) # abundance_table <- CompHeatmaps::create_abundance_table(processed_data) # CompHeatmaps::create_heatmap(abundance_table) # ``

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
