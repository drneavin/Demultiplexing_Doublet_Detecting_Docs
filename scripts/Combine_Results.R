#!/usr/bin/env Rscript

library(data.table)





# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, help="The output directory where results will be saved")
parser$add_argument("-d", "--demuxlet", required = FALSE, type = "character", default = NULL, help = "Path to demuxlet results. Only use this option if you want to include the demuxlet results.")
parser$add_argument("-f", "--freemuxlet", required = FALSE, type = "character", default=NULL, help = "Path to freemuxlet results. Only use this option if you want to include the freemuxlet results.")
parser$add_argument("-s", "--scSplit", required = FALSE, type="character", default=NULL, help="Path to scSplit results. Only use this option if you want to include the scSplit results.")
parser$add_argument("-u", "--souporcell", required = FALSE, type = "character", default=NULL, help = "Path to souporcell results. Only use this option if you want to include the souporcell results.")
parser$add_argument("-v", "--vireo", required = FALSE, type = "character", default=-NULL, help = "Path to vireo results. Only use this option if you want to include the vireo results.")
parser$add_argument("-d", "--DoubletDecon", required = FALSE, type = "character", default=NULL, help = "Path to DoubletDecon results. Only use this option if you want to include the DoubletDecon results.")
parser$add_argument("-c", "--DoubletDetection", required = FALSE, type = "character", default=NULL, help = "Path to DoubletDetection results. Only use this option if you want to include the DoubletDetection results.")
parser$add_argument("-f", "--DoubletFinder", required = FALSE, type = "character", default=NULL, help = "Path to DoubletFinder results. Only use this option if you want to include the DoubletFinder results.")
parser$add_argument("-i", "--scDblFinder", required = FALSE, type = "character", default=NULL, help = "Path to scDblFinder results. Only use this option if you want to include the scDblFinder results.")
parser$add_argument("-c", "--scds", required = FALSE, type = "character", default=NULL, help = "Path to scds results. Only use this option if you want to include the scds results.")
parser$add_argument("-r", "--scrublet", required = FALSE, type = "character", default=NULL, help = "Path to scrublet results. Only use this option if you want to include the scrublet results.")
parser$add_argument("-o", "--solo", required = FALSE, type = "character", default=NULL, help = "Path to solo results. Only use this option if you want to include the solo results.")
parser$add_argument("-m", "--method", required = FALSE, type = "character", default=NULL, help = "Combination method. Leave blank if you just want all the combine the methods into a single dataframe.")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()


### Check to make sure user has provided at least two softwares to combine ###
if (length(which(c(!is.null(args$demuxlet), !is.null(args$freemuxlet), !is.null(args$scSplit), !is.null(args$souporcell), !is.null(args$vireo),!is.null(args$DoubletDecon), !is.null(args$DoubletDetection), !is.null(args$DoubletFinder), !is.null(args$scDblFinder), !is.null(args$scds), !is.null(args$scrublet), !is.null(args$solo)))) < 2){
	message("You didn't provide at least two files to combine. Please try again - this time providing at least two inputs."
	q()
} else{
	### Make a list to store the results that will be combined in ###
	results_list <- list()


	if (!is.null(args$demuxlet)){
		message("Reading in demuxlet results.")
		results_list[["demuxlet"]] <- fread(args$demuxlet, delim = "\t")
	}

	if (!is.null(args$freemuxlet)){
		message("Reading in demuxlet results.")
		results_list[["freemuxlet"]] <- fread(args$freemuxlet, delim = "\t")
	}

	if (!is.null(args$scSplit)){
		message("Reading in demuxlet results.")
		results_list[["scSplit"]] <- fread(args$scSplit, delim = "\t")
	}

	if (!is.null(args$souporcell)){
		message("Reading in demuxlet results.")
		results_list[["souporcell"]] <- fread(args$souporcell, delim = "\t")
	}

	if (!is.null(args$vireo)){
		message("Reading in demuxlet results.")
		results_list[["vireo"]] <- fread(args$vireo, delim = "\t")
	}

	if (!is.null(args$DoubletDecon)){
		message("Reading in demuxlet results.")
		results_list[["DoubletDecon"]] <- fread(args$DoubletDecon, delim = "\t")
	}

	if (!is.null(args$DoubletDetection)){
		message("Reading in demuxlet results.")
		results_list[["DoubletDetection"]] <- fread(args$DoubletDetection, delim = "\t")
	}

	if (!is.null(args$DoubletFinder)){
		message("Reading in demuxlet results.")
		results_list[["DoubletFinder"]] <- fread(args$DoubletFinder, delim = "\t")
	}

	if (!is.null(args$scDblFinder)){
		message("Reading in demuxlet results.")
		results_list[["scDblFinder"]] <- fread(args$scDblFinder, delim = "\t")
	}

	if (!is.null(args$scds)){
		message("Reading in demuxlet results.")
		results_list[["scds"]] <- fread(args$scds, delim = "\t")
	}

	if (!is.null(args$scrublet)){
		message("Reading in demuxlet results.")
		results_list[["scrublet"]] <- fread(args$scrublet, delim = "\t")
	}

	if (!is.null(args$solo)){
		message("Reading in demuxlet results.")
		results_list[["solo"]] <- fread(args$solo, delim = "\t")
	}

}


#### Demuxlet
## Barcode = "BARCODE"
## singlet_doublet_droplet = "DROPLET.TYPE"
## assignment = "BEST.GUESS"
##	ie: 9_39,465_466,0.50
##  	9_39,9_39,0.00

#### Freemuxlet

### If not assigned individuals
### NEED to check!!!
## Barcode = "BARCODE"
## singlet_doublet_droplet = "DROPLET.TYPE"
## assignment = "BEST.GUESS"
##	ie: 12,1
##  	12,12

### If assigned individuals


#### scSplit

### If not assigned individuals
### NEED to check!!!
## Barcode = "Barcode"
## singlet_doublet_droplet in "Cluster"
## assignment in "Cluster"
##	ie: 12,1
##  	12,12

