#!/usr/bin/env Rscript

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, help="The output directory where results will be saved")
parser$add_argument("-t", "--tenX_matrix", required = TRUE, type = "character", help = "Path to the 10x filtered matrix directory or h5 file.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()


suppressMessages(suppressWarnings(library(scDblFinder)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(tidyverse)))



### Read in data as an sce object ###
if (file.exists(args$tenX_matrix)){
    print(paste0("Using the following counts: ", args$tenX_matrix))
    if endsWith(args$tenX_matrix, ".h5"){
        counts <- Read10X_h5(args$tenX_matrix, gene.column = 1)
    } else {
        counts <- Read10X(args$tenX_matrix, gene.column = 1)
    }
} else {
    print(paste0("Cannot find the counts matrix ", args$tenX_matrix))
}
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

