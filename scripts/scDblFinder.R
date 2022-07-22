#!/usr/bin/env Rscript

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))

# create parser object
parser <- ArgumentParser() 

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, help="The output directory where results will be saved")
parser$add_argument("-t", "--tenX_matrix", required = TRUE, type = "character", help = "Path to the 10x filtered matrix directory or h5 file.")
parser$add_argument("-b", "--barcodes_filtered", required = FALSE, type = "character", help = "Path to a list of filtered barcodes to use for doublet detection.")

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
    if (endsWith(args$tenX_matrix, ".h5")){
        counts <- Read10X_h5(args$tenX_matrix)
    } else {
        counts <- Read10X(args$tenX_matrix, gene.column = 1)
    }
} else {
    print(paste0("Cannot find the counts matrix ", args$tenX_matrix))
}


if (!is.null(args$barcodes_filtered)){
    if (file.exists(args$barcodes_filtered)){
        message("Reading in the filtered barcode list.")
        filtered_barcodes <- read_delim(args$barcodes_filtered, delim = "\t", col_names = "Barcodes")
        number_original_barcodes <- ncol(counts)
        
        message("\nFiltering counts for barcodes included in the filtered barcode list.\n")
        ### Filter for the barcodes list of interest ###
        if (is.list(counts)){
            barcodes_head <- head(colnames(counts[[grep("Gene", names(counts))]]))
            counts <- counts[[grep("Gene", names(counts))]][, colnames(counts[[grep("Gene", names(counts))]]) %in% filtered_barcodes$Barcodes]
        } else {
            barcodes_head <- head(colnames(counts))
            counts <- counts[, colnames(counts) %in% filtered_barcodes$Barcodes]
        }

        ### Provide user informatin with the number of original barcodes, those in the filter list and those after filtering
        number_filtered_barcodes <- ncol(counts)

        if (number_filtered_barcodes == 0){
            message("\nThere were no barcodes left in counts matrix after filtering on the barcode list.\nHere is a list of the first few barcodes in the filter barcode list provided:\n")
            message(paste(head(filtered_barcodes$Barcodes), colapse = "\n"))
            message("And here is a list of the first few barcodes in the original counts matrix:\n")
            message(paste(barcodes_head, colapse = "\n"))
            message("\nPlease check for any inconsistencies between the barcode formatting.\n\nExiting\n")
            q()
        } else {
            message(paste0("\nThe original number of barcodes in the counts matrix: ", number_original_barcodes, "\nThe number of barcodes in the user-provided barcode list :", nrow(filtered_barcodes), "\nThe number of barcodes after filtering for user-provided barcodes: ", number_filtered_barcodes, "\n"))
        }

    } else {
        message(paste0("Cannot find the provided filtered barcode list at ", args$barcodes_filtered, ".\n\nExiting."))
        q()
    }
}



if (is.list(counts)){
	sce <- SingleCellExperiment(list(counts=counts[[grep("Gene", names(counts))]]))
} else {
	sce <- SingleCellExperiment(list(counts=counts))
}

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

