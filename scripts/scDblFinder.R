#!/usr/bin/env Rscript

.libPaths("/usr/local/lib/R/site-library")
library(argparse)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, help="The output directory where results will be saved")
parser$add_argument("-t", "--tenX_matrix", required = TRUE, type = "character", help = "Path to the 10x filtered matrix directory.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()


library(scDblFinder)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)


print(paste0("Using the following counts directory: ", args$tenX_matrix))



### Read in data as an sce object ###
counts <- Read10X(args$tenX_matrix, gene.column = 1)
sce <- SingleCellExperiment(list(counts=counts))


## Calculate doublet ratio ###
doublet_ratio <- ncol(sce)/1000*0.008


### Calculate Singlets and Doublets ###
sce <- scDblFinder(sce, dbr=doublet_ratio)


 
### Make a dataframe of the results ###
results <- data.frame("Barcode" = rownames(colData(sce)), "scDblFinder_DropletType" = sce$scDblFinder.class, "scDblFinder_Score" = sce$scDblFinder.score)


write_delim(results, file = paste0(args$out,"/scDblFinder_doublets_singlets.tsv"), delim = "\t")

### Calculate number of doublets and singlets ###
summary <- as.data.frame(table(results$scDblFinder_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
write_delim(summary, paste0(args$out,"/scDblFinder_doublet_summary.tsv"), "\t")

