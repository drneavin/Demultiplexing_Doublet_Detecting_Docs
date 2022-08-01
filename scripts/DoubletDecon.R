#!/usr/bin/env Rscript

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, help="The output directory where results will be saved")
parser$add_argument("-s", "--seurat_object", required = TRUE, type = "character", help = "A QC, normalized seurat object with classifications/clusters as Idents() saved as an rds object.")
parser$add_argument("-g", "--num_genes", required = FALSE, type = "integer", default=50, help = "Number  of genes to use in \'Improved_Seurat_Pre_Process\' function.")
parser$add_argument("-r", "--rhop", required = FALSE, type="double", default=0.9, help="rhop to use in DoubletDecon - the number of SD from the mean to identify upper limit to blacklist")
parser$add_argument("-p", "--species", required = FALSE, type = "character", default="hsa", help = "The species of your sample. Can be scientific species name, KEGG ID, three letter species abbreviation, or NCBI ID.")
parser$add_argument("-n", "--nCores", required = FALSE, type = "double", default=-1, help = "The number of unique cores you would like to use to run DoubletDecon. By default, uses one less than available detected.")
parser$add_argument("-c", "--removeCC", required = FALSE, type = "logical", default=FALSE, help = "Whether to remove clusters enriched in cell cycle genes.")
parser$add_argument("-m", "--pmf", required = FALSE, type = "logical", default=TRUE, help = "Whether to use unique gene expression in doublet determination.")
parser$add_argument("-f", "--heatmap", required = FALSE, type = "logical", default=FALSE, help = "Whether to generate heatmaps.")
parser$add_argument("-t", "--centroids", required = FALSE, type = "logical", default=FALSE, help = "Whether to use centroids instead of medoids for doublet detecting.")
parser$add_argument("-d", "--num_doubs", required = FALSE, type = "integer", default=100, help = "The number of doublets to simulate for each cluster pair.")
parser$add_argument("-5", "--only50", required = FALSE, type = "logical", default=FALSE, help = "Whether to only compute doublets as 50:50 ratio. Default is to use other ratios as well.")
parser$add_argument("-u", "--min_uniq", required = FALSE, type = "integer", default=4, help = "Minimum number of unique genes to rescue a cluster identified as doublets.")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

suppressMessages(suppressWarnings(library(DoubletDecon)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(data.table)))


## make sure the directory exists ###
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)

## Read in Data ##
seurat <- readRDS(args$seurat_object)

## Preprocess ##
processed <- Improved_Seurat_Pre_Process(seurat, num_genes=args$num_genes, write_files=FALSE)

## Run Doublet Decon ##
results <- Main_Doublet_Decon(rawDataFile = processed$newExpressionFile, 
  groupsFile = processed$newGroupsFile, 
  filename = "DoubletDecon_results",
  location = paste0(args$out, "/"),
  fullDataFile = NULL, 
  removeCC = args$removeCC, 
  species = args$species, 
  rhop = args$rhop,
  write = TRUE, 
  PMF = args$pmf, 
  useFull = FALSE, 
  heatmap = args$heatmap, 
  centroids=args$centroids, 
  num_doubs=args$num_doubs, 
  only50=args$only50, 
  min_uniq=args$min_uniq, 
  nCores = args$nCores)




doublets <- read.table(paste0(args$out, "/Final_doublets_groups_DoubletDecon_results.txt"))
doublets$Barcode <- gsub("\\.", "-",rownames(doublets))
doublets$DoubletDecon_DropletType <- "doublet"
doublets$V1 <- NULL
doublets$V2 <- NULL


singlets <- read.table(paste0(args$out, "/Final_nondoublets_groups_DoubletDecon_results.txt"))
singlets$Barcode <- gsub("\\.", "-",rownames(singlets))
singlets$DoubletDecon_DropletType <- "singlet"
singlets$V1 <- NULL
singlets$V2 <- NULL

doublets_singlets <- rbind(singlets,doublets)

fwrite(doublets_singlets, paste0(args$out, "/DoubletDecon_doublets_singlets.tsv"), sep = "\t", append = FALSE)


### Make a summaruy of the number of singlets and doublets
summary <- as.data.frame(table(doublets_singlets$DoubletDecon_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
fwrite(summary, paste0(args$out,"/DoubletDecon_doublet_summary.tsv"), sep = "\t", append = FALSE)