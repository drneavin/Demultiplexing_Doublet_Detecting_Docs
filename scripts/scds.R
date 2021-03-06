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


suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(scds)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))


## Read in data
if (file.exists(args$tenX_matrix)){
    print(paste0("Using the following counts: ", args$tenX_matrix))
    if (endsWith(args$tenX_matrix, ".h5")){
        counts <- Read10X_h5(args$tenX_matrix, gene.column = 1)
    } else {
        counts <- Read10X(args$tenX_matrix, gene.column = 1)
    }
} else {
    print(paste0("Cannot find the counts matrix ", args$tenX_matrix))
}

if (is.list(counts)){
	sce <- SingleCellExperiment(list(counts=counts[[grep("Gene", names(counts))]]))
} else {
	sce <- SingleCellExperiment(list(counts=counts))
}

## Annotate doublet using binary classification based doublet scoring:
sce = bcds(sce, retRes = TRUE, estNdbl=TRUE)

## Annotate doublet using co-expression based doublet scoring:
try({
    sce = cxds(sce, retRes = TRUE, estNdbl=TRUE)
})

### If cxds worked, run hybrid, otherwise use bcds annotations
if ("cxds_score" %in% colnames(colData(sce))) {
    ## Combine both annotations into a hybrid annotation
    sce = cxds_bcds_hybrid(sce, estNdbl=TRUE)
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$hybrid_score, colData(sce)$hybrid_call))
} else {
    print("this pool failed cxds so results are just the bcds calls")
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$bcds_score, colData(sce)$bcds_call))
}

## Doublet scores are now available via colData:
colnames(Doublets) <- c("Barcode","scds_score","scds_DropletType")
Doublets$scds_DropletType <- gsub("FALSE","singlet",Doublets$scds_DropletType) 
Doublets$scds_DropletType <- gsub("TRUE","doublet",Doublets$scds_DropletType)

message("writing output")
write_delim(Doublets, paste0(args$out,"/scds_doublets_singlets.tsv"), "\t")


summary <- as.data.frame(table(Doublets$scds_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
write_delim(summary, paste0(args$out,"/scds_doublet_summary.tsv"), "\t")


