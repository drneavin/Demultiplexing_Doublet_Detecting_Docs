#!/usr/bin/env Rscript

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--assign", required = TRUE, help="The dropulation assignment results.")
parser$add_argument("-s", "--doublet", required = TRUE, type = "character", help = "The dropulation doublet calls.")
parser$add_argument("-c", "--out", required = TRUE, type = "logical", help = "Whether sctransform was used for normalization.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(data.table)))


createDonorMap<-function (args$assign, doubletFile, outMapFile, fdrThreshold=0.05, doubletPvalueThreshold=0.9) {
    a=read.table(args$assign, header=T, stringsAsFactors = F, sep="\t")
    b=read.table(doubletFile, header=T, stringsAsFactors = F, sep="\t")
    
    confidentAssignmentCells=a[a$FDR_pvalue<=fdrThreshold,"cell"]
    singletCells=b[b$doublet_pval<doubletPvalueThreshold, "cell"]
    confidentAssignmentSingletCells=intersect(confidentAssignmentCells, singletCells)
    doublets = a$cell[!(a$cell %in% confidentAssignmentSingletCells)]

    mapa=a[match(confidentAssignmentSingletCells, a$cell), c("cell", "bestLikelihood", "bestSample")]
    colnames(mapa) = c("Barcode","dropulation_Likelihood", "dropulation_Assignment")
    mapa$dropulation_DropletType = "singlet"

    mapb= b[match(doublets, b$cell), c("cell", "mixedSampleLikelihood", "mixedSample")]
    colnames(mapb) = c("Barcode", "dropulation_Likelihood","dropulation_Assignment")
    mapb$dropulation_Assignment = "doublet"
    mapb$dropulation_DropletType = "doublet"

    map = rbind(mapa,mapb)

    map = left_join(map, a[,c("cell", "num_snps", "num_umis")], by = c("Barcode" = "cell"))
    colnames(map) = gsub("num_snps", "dropulation_Nsnps", colnames(map)) %>%
                        gsub("num_umis", "dropulation_Numis", .)


    fwrite(map, outMapFile, row.names=F, col.names = T, quote=F, sep="\t")
}


createDonorMap(args$assign, args$doublet, parser$out)

