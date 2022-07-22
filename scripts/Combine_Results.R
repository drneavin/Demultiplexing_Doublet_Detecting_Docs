#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(future.apply)))
suppressMessages(suppressWarnings(library(ComplexUpset)))
suppressMessages(suppressWarnings(library(RColorBrewer)))


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, type = "character", help="The file where results will be saved")
parser$add_argument("-d", "--demuxlet", required = FALSE, type = "character", default = NULL, help = "Path to demuxlet results. Only use this option if you want to include the demuxlet results.")
parser$add_argument("-f", "--freemuxlet", required = FALSE, type = "character", default=NULL, help = "Path to freemuxlet results. Only use this option if you want to include the freemuxlet results.")
parser$add_argument("-g", "--freemuxlet_assignments", required = FALSE, type = "character", default=NULL, help = "Path to freemuxlet cluster-to-individual assignments. Only use this option if have used reference SNP genotypes to assign individuals to clusters for the freemuxlet results.")
parser$add_argument("-a", "--freemuxlet_correlation_limit", required = FALSE, type = "double", default=0.7, help = "The minimum correlation between the cluster and the individual SNP genotypes which should be considered as a valid assignment. If you want no limit, use 0. Default is 0.7.")
parser$add_argument("-s", "--scSplit", required = FALSE, type="character", default=NULL, help="Path to scSplit results. Only use this option if you want to include the scSplit results.")
parser$add_argument("-w", "--scSplit_assignments", required = FALSE, type="character", default=NULL, help="Path to scSplit cluster-to-individual assignments. Only use this option if you have used reference SNP genotypes to assign individuals to clusters for the scSplit results.")
parser$add_argument("-j", "--scSplit_correlation_limit", required = FALSE, type = "double", default=0.7, help = "The minimum correlation between the cluster and the individual SNP genotypes which should be considered as a valid assignment. If you want no limit, use 0. Default is 0.7.")
parser$add_argument("-u", "--souporcell", required = FALSE, type = "character", default=NULL, help = "Path to souporcell results. Only use this option if you want to include the souporcell results.")
parser$add_argument("-x", "--souporcell_assignments", required = FALSE, type = "character", default=NULL, help = "Path to souporcell cluster-to-individual assignments. Only use this option if you have used reference SNP genotypes to assign individuals to clusters for the souporcell results.")
parser$add_argument("-k", "--souporcell_correlation_limit", required = FALSE, type = "double", default=0.7, help = "The minimum correlation between the cluster and the individual SNP genotypes which should be considered as a valid assignment. If you want no limit, use 0. Default is 0.7.")
parser$add_argument("-v", "--vireo", required = FALSE, type = "character", default=NULL, help = "Path to vireo results. Only use this option if you want to include the vireo results.")
parser$add_argument("-e", "--DoubletDecon", required = FALSE, type = "character", default=NULL, help = "Path to DoubletDecon results. Only use this option if you want to include the DoubletDecon results.")
parser$add_argument("-t", "--DoubletDetection", required = FALSE, type = "character", default=NULL, help = "Path to DoubletDetection results. Only use this option if you want to include the DoubletDetection results.")
parser$add_argument("-i", "--DoubletFinder", required = FALSE, type = "character", default=NULL, help = "Path to DoubletFinder results. Only use this option if you want to include the DoubletFinder results.")
parser$add_argument("-n", "--scDblFinder", required = FALSE, type = "character", default=NULL, help = "Path to scDblFinder results. Only use this option if you want to include the scDblFinder results.")
parser$add_argument("-c", "--scds", required = FALSE, type = "character", default=NULL, help = "Path to scds results. Only use this option if you want to include the scds results.")
parser$add_argument("-r", "--scrublet", required = FALSE, type = "character", default=NULL, help = "Path to scrublet results. Only use this option if you want to include the scrublet results.")
parser$add_argument("-l", "--solo", required = FALSE, type = "character", default=NULL, help = "Path to solo results. Only use this option if you want to include the solo results.")
parser$add_argument("-b", "--ref", required = FALSE, type = "character", default=NULL, help = "Which demultiplexing software to use as a reference for individuals when you do not have assignment key for all demultiplexing method. Options are 'Demuxlet', 'Freemuxlet', 'scSplit', 'Souporcell' and 'Vireo'. If blank when assignment keys are missing, default softwares to use if present are Vireo, then Demuxlet, then Freemuxlet, then Souporcell, then scSplit.")
parser$add_argument("-p", "--pct_agreement", required = FALSE, type = "double", default=0.7, help = "The proportion of a cluster that match the 'ref' assignment to assign that cluster the individual assignment from the reference. Can be between 0.5 and 1. Default is 0.9.")
parser$add_argument("-m", "--method", required = FALSE, type = "character", default=NULL, help = "Combination method. Options are 'MajoritySinglet'. 'AtLeastHalfSinglet', 'AnySinglet' or 'AnyDoublet'. We have found that 'MajoritySinglet' provides the most accurate results in most situations and therefore recommend this method. See https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/CombineResults.html for detailed explanation of each intersectional method. Leave blank if you just want all the softwares to be merged into a single dataframe.")
                                        

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()


### Functions ###
which.max.simple=function(x,na.rm=TRUE,tie_value="NA"){
	if(na.rm){
		x=x[!is.na(x)]
	}
	if(length(x)==0){
		return(NA)
	}
	maxval=max(x)
	if(is.na(maxval)){
		return(NA)
	}
	if(sum(x %in% maxval) > 1){
		# Ties exist, figure out what to do with them.
		if(tie_value=="NA"){
			return(NA)
		}
		if(tie_value=="random"){
			tie_postions=which(x==maxval)
			return(sample(tie_postions,size=1))
		}
		if(tie_value=="first"){
			tie_postions=which(x==maxval)
			return(tie_postions[1])
		}
	} else{
		return(which.max(x))
	}
}


### Check to be sure that output directory already exists or create it for the user ###
if (!dir.exists(gsub(basename(args$out), "", args$out))){
	message(paste0("Creating directory as it doesn't exist. Complete path: ", gsub(basename(args$out), "", args$out)))
	dir.create(gsub(basename(args$out), "", args$out), recursive = TRUE)
}


### Check to make sure user has provided at least two softwares to combine ###
if (length(which(c(!is.null(args$demuxlet), !is.null(args$freemuxlet), !is.null(args$scSplit), !is.null(args$souporcell), !is.null(args$vireo),!is.null(args$DoubletDecon), !is.null(args$DoubletDetection), !is.null(args$DoubletFinder), !is.null(args$scDblFinder), !is.null(args$scds), !is.null(args$scrublet), !is.null(args$solo)))) < 1){
	message("You didn't provide any software results to combine. Please try again - this time providing at least one inputs.")
	q()
} else{
	### Make a list to store the results that will be combined in ###
	results_list <- list()

	if (!is.null(args$demuxlet)){

		### Read in Demuxlet Results ###
		message("Reading in demuxlet results.")
		if (file_test("-f", args$demuxlet)){
				results_list[["demuxlet"]] <- fread(args$demuxlet, sep = "\t")
		} else {
			demuxlet <- list.files(args$demuxlet, pattern = ".best")
			if (length(demuxlet) == 1){
				results_list[["demuxlet"]] <- fread(paste0(args$demuxlet, "/", demuxlet), sep = "\t")
			} else {
				message("Can't read in the demuxlet report either from the file you provided or to find the '*.best' file in the directory you provided. Please double check your path and provide the full path to the demuxlet file. Exiting.")
				q()
			}
		}

		### Update dataframe to just be categories of interest
		results_list[["demuxlet"]] <- results_list[["demuxlet"]][,c("BARCODE", "DROPLET.TYPE", "BEST.GUESS")]
		colnames(results_list[["demuxlet"]]) <- c("Barcode", "Demuxlet_DropletType", "Demuxlet_Individual_Assignment")

		### update Demuxlet_DropletType  column SNG and DBL ###
		results_list[["demuxlet"]]$Demuxlet_DropletType  <- gsub("SNG", "singlet", results_list[["demuxlet"]]$Demuxlet_DropletType) %>% gsub("DBL", "doublet", .) %>% gsub("AMB", "unassigned", .)

		### Update Individual column to be individual or doublet or ambiguous/unassigned ###
		results_list[["demuxlet"]]$Demuxlet_Individual_Assignment <- ifelse(results_list[["demuxlet"]]$Demuxlet_DropletType  == "doublet", "doublet",
																		ifelse(results_list[["demuxlet"]]$Demuxlet_DropletType  == "unassigned", "unassigned", gsub(",.+,.\\..+","",results_list[["demuxlet"]]$Demuxlet_Individual_Assignment)))

	}


	if (!is.null(args$freemuxlet)){
		message("Reading in freemuxlet results.")
		if (file_test("-f", args$freemuxlet)){
			results_list[["Freemuxlet"]] <- fread(args$freemuxlet, sep = "\t")
		} else {
			freemuxlet <- list.files(args$freemuxlet, pattern = ".samples.gz")
			if (length(freemuxlet) == 1){
					results_list[["Freemuxlet"]] <- fread(paste0(args$freemuxlet, "/", freemuxlet), sep = "\t")
			} else {
				message("Can't read in the freemuxlet report either from the file you provided or to find the '*.samples.gz' file in the directory you provided. Please double check your path and provide the full path to the freemuxlet file. Exiting.")
				q()
			}
		}

		### Update dataframe to just be categories of interest
		results_list[["Freemuxlet"]] <- results_list[["Freemuxlet"]][,c("BARCODE", "DROPLET.TYPE", "BEST.GUESS")]
		colnames(results_list[["Freemuxlet"]]) <- c("Barcode", "Freemuxlet_DropletType", "Freemuxlet_Cluster")

		### update Demuxlet_DropletType  column SNG and DBL ###
		results_list[["Freemuxlet"]]$Freemuxlet_DropletType  <- gsub("SNG", "singlet", results_list[["Freemuxlet"]]$Freemuxlet_DropletType ) %>% gsub("DBL", "doublet", .) %>% gsub("AMB", "unassigned", .)

		### Update Individual column to be individual or doublet or ambiguous/unassigned ###
		results_list[["Freemuxlet"]]$Freemuxlet_Cluster <- ifelse(results_list[["Freemuxlet"]]$Freemuxlet_DropletType  == "doublet", "doublet",
																		ifelse(results_list[["Freemuxlet"]]$Freemuxlet_DropletType  == "unassigned", "unassigned", gsub(",.+","",results_list[["Freemuxlet"]]$Freemuxlet_Cluster)))

	}

	if (!is.null(args$scSplit)){
		message("Reading in scSplit results.")
		if (file_test("-f", args$scSplit)){
			results_list[["scSplit"]] <- fread(args$scSplit, sep = "\t")
		} else {
			tryCatch({
				results_list[["scSplit"]] <- fread(paste0(args$scSplit, "/scSplit_result.csv"), sep = "\t")
			},
			error = function(e) {
				message(e)
				message("Can't read in the scSplit report either from the file you provided or to find the 'scSplit_result.csv' file in the directory you provided. Please double check your path and provide the full path to the scSplit file. Exiting.")
				q()
			})
		}

		### Update dataframe to just be categories of interest
		results_list[["scSplit"]][, c("scSplit_DropletType", "scSplit_Cluster") := tstrsplit(Cluster, "-", fixed=TRUE)]
		results_list[["scSplit"]]$Cluster <- NULL
		results_list[["scSplit"]]$scSplit_DropletType  <- gsub("SNG", "singlet", results_list[["scSplit"]]$scSplit_DropletType) %>%
															gsub("DBL", "doublet", .) 
		results_list[["scSplit"]]$scSplit_Cluster <- ifelse(results_list[["scSplit"]]$scSplit_DropletType  == "doublet", "doublet", results_list[["scSplit"]]$scSplit_Cluster)

	}

	if (!is.null(args$souporcell)){
		message("Reading in souporcell results.")
		if (file_test("-f", args$souporcell)){
			results_list[["Souporcell"]] <- fread(args$souporcell, sep = "\t")
		} else {
			tryCatch({
				results_list[["Souporcell"]] <- fread(paste0(args$souporcell, "/clusters.tsv"), sep = "\t")
			},
			error = function(e) {
				message(e)
				message("Can't read in the souporcell report either from the file you provided or to find the 'clusters.tsv' file in the directory you provided. Please double check your path and provide the full path to the souporcell file. Exiting.")
				q()
			})
		}

		### Update dataframe to just be categories of interest
		results_list[["Souporcell"]] <- results_list[["Souporcell"]][,c("barcode", "status","assignment")]
		colnames(results_list[["Souporcell"]]) <- c("Barcode", "Souporcell_DropletType", "Souporcell_Cluster")
		results_list[["Souporcell"]]$Souporcell_Cluster <- ifelse(results_list[["Souporcell"]]$Souporcell_DropletType  == "doublet", "doublet", ifelse(results_list[["Souporcell"]]$Souporcell_DropletType  == "unassigned", "unassigned", results_list[["Souporcell"]]$Souporcell_Cluster))

	}

	if (!is.null(args$vireo)){
		message("Reading in vireo results.")
		if (file_test("-f", args$vireo)){
			results_list[["Vireo"]] <- fread(args$vireo, sep = "\t")
		} else {
			tryCatch({
				results_list[["Vireo"]] <- fread(paste0(args$vireo, "/donor_ids.tsv"), sep = "\t")
			},
			error = function(e) {
				message(e)
				message("Can't read in the vireo report either from the file you provided or to find the 'donor_ids.tsv' file in the directory you provided. Please double check your path and provide the full path to the vireo file. Exiting.")
				q()
			})
		}
		
		### Update dataframe to just be categories of interest
		results_list[["Vireo"]] <- results_list[["Vireo"]][,c("cell", "donor_id")]
		colnames(results_list[["Vireo"]]) <- c("Barcode", "Vireo_Individual_Assignment")
		results_list[["Vireo"]]$Vireo_DropletType  <- ifelse(results_list[["Vireo"]]$Vireo_Individual_Assignment == "doublet", "doublet", ifelse(results_list[["Vireo"]]$Vireo_Individual_Assignment == "unassigned", "unassigned", "singlet"))

	}

	if (!is.null(args$DoubletDecon)){
		message("Reading in DoubletDecon results.")
		if (file_test("-f", args$DoubletDecon)){
			results_list[["DoubletDecon"]] <- fread(args$DoubletDecon, sep = "\t")
		} else {
			tryCatch({
				results_list[["DoubletDecon"]] <- fread(paste0(args$DoubletDecon, "/DoubletDecon_doublets_singlets.tsv"), sep = "\t")
			},
			error = function(e) {
				message(e)
				message("Can't read in the DoubletDecon report either from the file you provided or to find the 'DoubletDecon_doublets_singlets.tsv' file in the directory you provided. Please double check your path and provide the full path to the DoubletDecon file. Exiting.")
				q()
			})
		}
	}

	if (!is.null(args$DoubletDetection)){
		message("Reading in DoubletDetection results.")
		if (file_test("-f", args$DoubletDetection)){
			results_list[["DoubletDetection"]] <- fread(args$DoubletDetection, sep = "\t")
		} else {
			tryCatch({
				results_list[["DoubletDetection"]] <- fread(paste0(args$DoubletDetection, "/DoubletDetection_doublets_singlets.tsv"), sep = "\t")
			},
			error = function(e) {
				message(e)
				message("Can't read in the DoubletDetection report either from the file you provided or to find the 'DoubletDetection_doublets_singlets.tsv' file in the directory you provided. Please double check your path and provide the full path to the DoubletDetection file. Exiting.")
				q()
			})
		}
	}

	if (!is.null(args$DoubletFinder)){
		message("Reading in DoubletFinder results.")
		if (file_test("-f", args$DoubletFinder)){
			results_list[["DoubletFinder"]] <- fread(args$DoubletFinder, sep = "\t")
		} else {
			tryCatch({
				results_list[["DoubletFinder"]] <- fread(paste0(args$DoubletFinder, "/DoubletFinder_doublets_singlets.tsv"), sep = "\t")
			},
			error = function(e) {
				message("Can't read in the DoubletFinder report either from the file you provided or to find the 'DoubletFinder_doublets_singlets.tsv' file in the directory you provided. Please double check your path and provide the full path to the DoubletFinder file. Exiting.")
				q()
			})
		}
	}

	if (!is.null(args$scDblFinder)){
		message("Reading in scDblFinder results.")
		if (file_test("-f", args$scDblFinder)){
			results_list[["scDblFinder"]] <- fread(args$scDblFinder, sep = "\t")
		} else {
			tryCatch({
				results_list[["scDblFinder"]] <- fread(paste0(args$scDblFinder, "/scDblFinder_doublets_singlets.tsv"), sep = "\t")
			},
			error = function(e) {
				message("Can't read in the scDblFinder report either from the file you provided or to find the 'scDblFinder_doublets_singlets.tsv' file in the directory you provided. Please double check your path and provide the full path to the scDblFinder file. Exiting.")
				q()
			})
		}
	}

	if (!is.null(args$scds)){
		message("Reading in scds results.")
		if (file_test("-f", args$scds)){
			results_list[["scds"]] <- fread(args$scds, sep = "\t")
		} else {
			tryCatch({
				results_list[["scds"]] <- fread(paste0(args$scds, "/scds_doublets_singlets.tsv"), sep = "\t")
			},
			error = function(e) {
				message("Can't read in the scds report either from the file you provided or to find the 'scds_doublets_singlets.tsv' file in the directory you provided. Please double check your path and provide the full path to the scds file. Exiting.")
				q()
			})
		}
	}

	if (!is.null(args$scrublet)){
		message("Reading in scrbulet results.")
		if (file_test("-f", args$scrublet)){
			results_list[["scrublet"]] <- fread(args$scrublet, sep = "\t")
		} else {
			tryCatch({
				results_list[["scrublet"]] <- fread(paste0(args$scrublet, "/scrublet_results.tsv"), sep = "\t")
			},
			error = function(e) {
				message("Can't read in the scrublet report either from the file you provided or to find the 'scrublet_results.tsv' file in the directory you provided. Please double check your path and provide the full path to the scrublet file. Exiting.")
				q()
			})
		}
	}

	if (!is.null(args$solo)){
		message("Reading in solo results.")
		if (file_test("-f", args$solo)){
			results_list[["solo"]] <- fread(args$solo, sep = "\t")
		} else {
			tryCatch({
				results_list[["solo"]] <- fread(paste0(args$solo, "/solo_results.tsv"), sep = "\t")
			},
			error = function(e) {
				message("Can't read in the solo report either from the file you provided or to find the 'solo_results.tsv' file in the directory you provided. Please double check your path and provide the full path to the solo file. Exiting.")
				q()
			})
		}
	}



	### Check for cluster-to-individual assignments + read in ###
	if (!is.null(args$freemuxlet_assignments) | !is.null(args$scSplit_assignments) | !is.null(args$souporcell_assignments)){
		### Make a list to store the results that will be combined in ###
		results_assignments_list <- list()

		if (!is.null(args$freemuxlet_assignments)){
			message("Reading in freemuxlet cluster-to-indiviudal assignments.")
			if (file_test("-f", args$freemuxlet_assignments)){
				results_assignments_list[["Freemuxlet"]] <- fread(args$freemuxlet_assignments, sep = "\t")
			} else {
				tryCatch({
					results_assignments_list[["Freemuxlet"]] <- fread(paste0(args$freemuxlet_assignments, "/Genotype_ID_key.txt"), sep = "\t")
				},
				error = function(e) {
					message("Can't read in the freemuxlet cluster-to-indiviudal assignments report either from the file you provided or to find the 'Genotype_ID_key.tsv' file in the directory you provided. Please double check your path and provide the full path to the freemuxlet cluster-to-indiviudal key file. Exiting.")
					q()
				})
			}


			### Check for column names ###
			if (all(!(grepl("Genotype_ID",colnames(results_assignments_list[["Freemuxlet"]]))) & !(grepl("Cluster_ID", colnames(results_assignments_list[["Freemuxlet"]]))))){
				message("Didn't find 'Genotype_ID' and 'Cluster_ID' in the columns of your freemuxlet cluster-to-indiviudal assignments file. Will use first column as the individual ID and the second column will be the cluster ID.")
			} else {
				if ("Correlation" %in% colnames(results_assignments_list[["Freemuxlet"]])){
					temp <- results_assignments_list[["Freemuxlet"]][,c("Genotype_ID", "Cluster_ID", "Correlation")]
					results_assignments_list[["Freemuxlet"]] <- temp
				}
			}
			colnames(results_assignments_list[["Freemuxlet"]])[1:2] <- c("Freemuxlet_Individual_Assignment", "Freemuxlet_Cluster")
			results_assignments_list[["Freemuxlet"]]$Freemuxlet_Cluster <- as.character(results_assignments_list[["Freemuxlet"]]$Freemuxlet_Cluster)

			### remove "CLUST" from assignments ###
			results_assignments_list[["Freemuxlet"]]$Freemuxlet_Cluster <- gsub("CLUST", "", results_assignments_list[["Freemuxlet"]]$Freemuxlet_Cluster)
			results_assignments_list[["Freemuxlet"]] <- results_assignments_list[["Freemuxlet"]][Freemuxlet_Cluster != "unassigned"]
			
			if ("Correlation" %in% colnames(results_assignments_list[["Freemuxlet"]])){
				results_assignments_list[["Freemuxlet"]] <- results_assignments_list[["Freemuxlet"]][Correlation > args$freemuxlet_correlation_limit]
				results_assignments_list[["Freemuxlet"]]$Correlation <- NULL
			}

			### Join
			message("Adding assignments to freeuxlet dataframe")
			results_list[["Freemuxlet"]] <- results_assignments_list[["Freemuxlet"]][results_list[["Freemuxlet"]], on="Freemuxlet_Cluster"]  

			results_list[["Freemuxlet"]]$Freemuxlet_Individual_Assignment <- as.character(ifelse(results_list[["Freemuxlet"]]$Freemuxlet_DropletType == "doublet", "doublet", ifelse(results_list[["Freemuxlet"]]$Freemuxlet_DropletType == "unassigned", "unassigned", results_list[["Freemuxlet"]]$Freemuxlet_Individual_Assignment)))
			results_list[["Freemuxlet"]]$Freemuxlet_Individual_Assignment <- as.character(ifelse(is.na(results_list[["Freemuxlet"]]$Freemuxlet_Individual_Assignment), results_list[["Freemuxlet"]]$Freemuxlet_Cluster, results_list[["Freemuxlet"]]$Freemuxlet_Individual_Assignment))


			### check for excessive NA after joining ###
			if (length(which(is.na(results_list[["Freemuxlet"]]$Freemuxlet_Individual_Assignment)))/nrow(results_list[["Freemuxlet"]]) >= 0.5){
				message("WARNING: There are more than 50% NA assignments after joining the cluster-to-individual assignments for freemuxlet. Please double check you files if this is unexpected.")
			}
		}


		if (!is.null(args$scSplit_assignments)){
			message("Reading in scSplit_assignments cluster-to-indiviudal assignments.")
			if (file_test("-f", args$scSplit_assignments)){
				results_assignments_list[["scSplit"]] <- fread(args$scSplit_assignments, sep = "\t")
			} else {
				tryCatch({
					results_assignments_list[["scSplit"]] <- fread(paste0(args$scSplit_assignments, "/Genotype_ID_key.txt"), sep = "\t")
				},
				error = function(e) {
					message("Can't read in the scSplit cluster-to-indiviudal assignments report either from the file you provided or to find the 'Genotype_ID_key.tsv' file in the directory you provided. Please double check your path and provide the full path to the scSplit cluster-to-indiviudal key file. Exiting.")
					q()
				})
			}


			### Check for column names ###
			if (all(!(grepl("Genotype_ID",colnames(results_assignments_list[["scSplit"]]))) & !(grepl("Cluster_ID", colnames(results_assignments_list[["scSplit"]]))))){
				message("Didn't find 'Genotype_ID' and 'Cluster_ID' in the columns of your scSplit cluster-to-indiviudal assignments. Will use first column as the individual ID and the second column will be the cluster ID.")
			} else {
				if ("Correlation" %in% colnames(results_assignments_list[["scSplit"]])){
					temp <- results_assignments_list[["scSplit"]][,c("Genotype_ID", "Cluster_ID", "Correlation")]
					results_assignments_list[["scSplit"]] <- temp
				}
					colnames(results_assignments_list[["scSplit"]])[1:2] <- c("scSplit_Individual_Assignment", "scSplit_Cluster")
			}
			colnames(results_assignments_list[["scSplit"]])[1:2] <- c("scSplit_Individual_Assignment", "scSplit_Cluster")
			results_assignments_list[["scSplit"]]$scSplit_Cluster <- as.character(results_assignments_list[["scSplit"]]$scSplit_Cluster)
			results_assignments_list[["scSplit"]] <- results_assignments_list[["scSplit"]][scSplit_Cluster != "unassigned"]

			if ("Correlation" %in% colnames(results_assignments_list[["scSplit"]])){
				results_assignments_list[["scSplit"]] <- results_assignments_list[["scSplit"]][Correlation > args$scSplit_correlation_limit]
				results_assignments_list[["scSplit"]]$Correlation <- NULL
			}

			### Join
			message("Adding assignments to scSplit dataframe")
			results_list[["scSplit"]] <- results_assignments_list[["scSplit"]][results_list[["scSplit"]], on="scSplit_Cluster"]  

			results_list[["scSplit"]]$scSplit_Individual_Assignment <- as.character(ifelse(results_list[["scSplit"]]$scSplit_DropletType == "doublet", "doublet", results_list[["scSplit"]]$scSplit_Individual_Assignment))
			results_list[["scSplit"]]$scSplit_Individual_Assignment <- as.character(ifelse(is.na(results_list[["scSplit"]]$scSplit_Individual_Assignment), results_list[["scSplit"]]$scSplit_Cluster, results_list[["scSplit"]]$scSplit_Individual_Assignment))

			### check for excessive NA after joining ###
			if (length(which(is.na(results_list[["scSplit"]]$Freemuxlet_Individual_Assignment)))/nrow(results_list[["scSplit"]]) >= 0.5){
				message("WARNING: There are more than 50% NA assignments after joining the cluster-to-individual assignments for scSplit. Please double check you files if this is unexpected.")
			}
		}

		if (!is.null(args$souporcell_assignments)){
			message("Reading in souporcell_assignments cluster-to-indiviudal assignments.")
			if (file_test("-f", args$souporcell_assignments)){
				results_assignments_list[["Souporcell"]] <- fread(args$souporcell_assignments, sep = "\t")
			} else {
				tryCatch({
					results_assignments_list[["Souporcell"]] <- fread(paste0(args$souporcell_assignments, "/Genotype_ID_key.txt"), sep = "\t")
				},
				error = function(e) {
					message("Can't read in the souporcell cluster-to-indiviudal assignments report either from the file you provided or to find the 'Genotype_ID_key.tsv' file in the directory you provided. Please double check your path and provide the full path to the souporcell cluster-to-indiviudal key file. Exiting.")
					q()
				})
			}


			### Check for column names ###
			if (all(!(grepl("Genotype_ID",colnames(results_assignments_list[["Souporcell"]]))) & !(grepl("Cluster_ID", colnames(results_assignments_list[["Souporcell"]]))))){
				message("Didn't find 'Genotype_ID' and 'Cluster_ID' in the columns of your souporcell cluster-to-indiviudal assignments. Will use first column as the individual ID and the second column will be the cluster ID.")
			} else {
				if ("Correlation" %in% colnames(results_assignments_list[["Souporcell"]])){
					temp <- results_assignments_list[["Souporcell"]][,c("Genotype_ID", "Cluster_ID", "Correlation")]
					results_assignments_list[["Souporcell"]] <- temp
				} 
					colnames(results_assignments_list[["Souporcell"]])[1:2] <- c("Souporcell_Individual_Assignment", "Souporcell_Cluster")
			}
			colnames(results_assignments_list[["Souporcell"]])[1:2] <- c("Souporcell_Individual_Assignment", "Souporcell_Cluster")
			results_assignments_list[["Souporcell"]]$Souporcell_Cluster <- as.character(results_assignments_list[["Souporcell"]]$Souporcell_Cluster)
			results_assignments_list[["Souporcell"]] <- results_assignments_list[["Souporcell"]][Souporcell_Cluster != "unassigned"]
			
			if ("Correlation" %in% colnames(results_assignments_list[["Souporcell"]])){
				results_assignments_list[["Souporcell"]] <- results_assignments_list[["Souporcell"]][Correlation > args$souporcell_correlation_limit]
				results_assignments_list[["Souporcell"]]$Correlation <- NULL
			}

			### Join
			message("Adding assignments to souporcell dataframe")
			results_list[["Souporcell"]] <- results_assignments_list[["Souporcell"]][results_list[["Souporcell"]], on="Souporcell_Cluster"]  

			results_list[["Souporcell"]]$Souporcell_Individual_Assignment <- as.character(ifelse(results_list[["Souporcell"]]$Souporcell_DropletType == "doublet", "doublet", ifelse(results_list[["Souporcell"]]$Souporcell_DropletType == "unassigned", "unassigned", results_list[["Souporcell"]]$Souporcell_Individual_Assignment)))
			results_list[["Souporcell"]]$Souporcell_Individual_Assignment <- as.character(ifelse(is.na(results_list[["Souporcell"]]$Souporcell_Individual_Assignment), results_list[["Souporcell"]]$Souporcell_Cluster, results_list[["Souporcell"]]$Souporcell_Individual_Assignment))

			### check for excessive NA after joining ###
			if (length(which(is.na(results_list[["Souporcell"]]$Freemuxlet_Individual_Assignment)))/nrow(results_list[["Souporcell"]]) >= 0.5){
				message("WARNING: There are more than 50% NA assignments after joining the cluster-to-individual assignments for souporcell. Please double check your files if this is unexpected.")
			}
			print(results_list[["Souporcell"]] )
		}

	}

	message("\nCombining results.\n")
	combined_results <- reduce(results_list, left_join)



	### Check for softwares with clusters but not assignments to standardize across softwares ###
	cluster_no_assign <- c()
	update_assignments <- FALSE


	if (length(c(args$demuxlet, args$freemuxlet, args$scSplit, args$souporcell, args$vireo)) > 1) {


		if (("Souporcell_Cluster" %in% colnames(combined_results)) & !("Souporcell_Individual_Assignment" %in% colnames(combined_results))){
			cluster_no_assign <- c(cluster_no_assign, "Souporcell")
			update_assignments <- TRUE
		} else if (!is.null(args$souporcell_assignments)){
			if ("Souporcell_Individual_Assignment" %in% colnames(combined_results)) {
				if (length(unique(results_assignments_list[["Souporcell"]]$Souporcell_Cluster)) < length(unique(combined_results$Souporcell_Individual_Assignment[!(combined_results$Souporcell_Individual_Assignment %in% c("doublet", "unassigned"))]))) {
					update_assignments <- TRUE
				}
			}
		}

		if (("scSplit_Cluster" %in% colnames(combined_results)) & !("scSplit_Individual_Assignment" %in% colnames(combined_results))){
			cluster_no_assign <- c(cluster_no_assign, "scSplit")
			update_assignments <- TRUE
		}else if (!is.null(args$scSplit_assignments)){
			if ("scSplit_Individual_Assignment" %in% colnames(combined_results)) {
				if (length(unique(results_assignments_list[["scSplit"]]$scSplit_Cluster)) < length(unique(combined_results$scSplit_Individual_Assignment[!(combined_results$scSplit_Individual_Assignment %in% c("doublet", "unassigned"))]))){
					update_assignments <- TRUE
				}
			}
		}	

		if (("Freemuxlet_Cluster" %in% colnames(combined_results)) & !("Freemuxlet_Individual_Assignment" %in% colnames(combined_results))){
			cluster_no_assign <- c(cluster_no_assign, "Freemuxlet")
			update_assignments <- TRUE
		} else if (!is.null(args$freemuxlet_assignments)){
			if ("Freemuxlet_Individual_Assignment" %in% colnames(combined_results)) {
				if (length(unique(results_assignments_list[["Freemuxlet"]]$Freemuxlet_Cluster)) < length(unique(combined_results$Freemuxlet_Individual_Assignment[!(combined_results$Freemuxlet_Individual_Assignment %in% c("doublet", "unassigned"))]))) {
					update_assignments <- TRUE
					cluster_no_assign <- c(cluster_no_assign, "Freemuxlet")
				}
			}
		}
	}


	if (update_assignments){

		message("\nConverting demultiplexing clusters without individual assignments to common ids.\n")


		cluster_assign_cols <- c(paste0(cluster_no_assign,"_Cluster"), grep("Individual_Assignment", colnames(combined_results), value = TRUE))
		if (!is.null(args$ref)){
			if (!(args$ref %in% colnames(combined_results))){
				ref <- args$ref
			} else {
				message("Couldn't find data for ", args$ref, " from your input. Please double check that you provided results from ", args$ref, " when you submitted this code.")
			}
		} else {
			if ("Vireo_Individual_Assignment" %in% colnames(combined_results)) {
				ref <- "Vireo"
			} else if ("Demuxlet_Individual_Assignment" %in% colnames(combined_results)) {
				ref <- "Demuxlet"
			} else if ("Freemuxlet_Individual_Assignment" %in% colnames(combined_results)) {
				ref <- "Freemuxlet"
			} else if ("Souporcell_Individual_Assignment" %in% colnames(combined_results)) {
				ref <- "Souporcell"
			} else if ("scSplit_Individual_Assignment" %in% colnames(combined_results)) {
				ref <- "scSplit"
			} else if ("Freemuxlet_Cluster" %in% colnames(combined_results)) {
				ref <- "Freemuxlet"
			} else if ("Souporcell_Cluster" %in% colnames(combined_results)) {
				ref <- "Souporcell"
			} else if ("scSplit_Cluster" %in% colnames(combined_results)) {
				ref <- "scSplit"
			}
			message("Using ", ref, " as individual assignment reference for assignment standardization.")
		}

		### Fix cluster_assign_cols if was one of the softwares without an individual assignment
		cluster_assign_cols <- cluster_assign_cols[!(cluster_assign_cols %in% paste0(ref,"_Cluster"))]

		### Generate IDs if scSplit, Souporcell or freemuxlet that just have weird numbers etc as clusters
		if ((paste0(ref, "_Cluster") %in% colnames(combined_results)) & !(paste0(ref,"_Individual_Assignment") %in% colnames(combined_results))){
			clusters <- data.table(Cluster = unique(combined_results[[paste0(ref,"_Cluster")]])[!unique(combined_results[[paste0(ref,"_Cluster")]]) %in% c("doublet", "unassigned")])
			colnames(clusters) <- paste0(ref, "_Cluster")
			clusters[[paste0(ref,"_Individual_Assignment")]] <- paste0("donor",seq(1:nrow(clusters)))
			combined_results <- clusters[combined_results, on = paste0(ref, "_Cluster")]
		}

		assign_table <- list()
		ids <- list()
		assign_table_top <- list()
		cluster_df <- list()
		combined_results_subset <- list()

		for (software in cluster_no_assign){
			### Get a vector of the unique individual assignments
			### first check if assigned some by correlations with genotypes (ie genotyped some in pool but not all) and subset out the already assigned
			if (paste0(software, "_Individual_Assignment") %in% colnames(combined_results)){
				combined_results_subset[[software]] <- combined_results[combined_results[[paste0(software, "_Individual_Assignment")]] %in% combined_results[[unique(paste0(software, "_Cluster"))]]]

				assign_table[[software]] <- combined_results_subset[[software]][,list(Total=.N) , c(paste0(software,"_Cluster"))][combined_results_subset[[software]][,list(Count=.N) ,c(paste0(software,"_Cluster"), paste0(ref,"_Individual_Assignment"))], on = paste0(software,"_Cluster")]

			} else {
				assign_table[[software]] <- combined_results[,list(Total=.N) , c(paste0(software,"_Cluster"))][combined_results[,list(Count=.N) ,c(paste0(software,"_Cluster"), paste0(ref,"_Individual_Assignment"))], on = paste0(software,"_Cluster")]
			}

			assign_table[[software]]$Proportion <- assign_table[[software]]$Count/assign_table[[software]]$Total
			assign_table[[software]]$Total <- NULL
			assign_table[[software]] <- assign_table[[software]][rev(order(Count))]


			### get just the top occurences of those assignments from the table
			### First message if any of the softwares will lose a cluster
			if (nrow(assign_table[[software]][Proportion >= args$pct_agreement]) == 0 | !all(assign_table[[software]][Proportion >= args$pct_agreement][[paste0(software,"_Cluster")]] %in% unique(assign_table[[software]][[paste0(software,"_Cluster")]]))){
				message(paste0("WARNING: The following clusters/individuals from ", software, " had a lower percent agreement than the threshold set when executing this script (", args$pct_agreement, "):\n", paste(colapse = "\n", unique(assign_table[[software]][[paste0(software,"_Cluster")]])[!(unique(assign_table[[software]][[paste0(software,"_Cluster")]]) %in% assign_table[[software]][Proportion >= args$pct_agreement][[paste0(software,"_Cluster")]])]), "\n\nPlease check whether this threshold is appropriate for this dataset and rerun with a different threshold for --pct_agreement if it seems too stringent."))
			}


			assign_table_top[[software]] <- assign_table[[software]][Proportion >= args$pct_agreement]
			cluster_df[[software]] <- assign_table_top[[software]][, .SD, .SD = c(paste0(software, "_Cluster"), paste0(ref, "_Individual_Assignment"))]
			colnames(cluster_df[[software]]) <- c(paste0(software, "_Cluster"), paste0(software, "_Individual_Assignment"))


			### Add results to the correct datatable if doesn't exist and update if already exists
			if (!nrow(assign_table[[software]][Proportion >= args$pct_agreement]) == 0){
				if (paste0(software, "_Individual_Assignment") %in% colnames(results_list[[software]])){
					for (row in 1:nrow(cluster_df[[software]])){
						results_list[[software]][[paste0(software,"_Individual_Assignment")]] <- gsub((paste0("^", cluster_df[[software]][row,][[paste0(software, "_Cluster")]], "$")), as.character(cluster_df[[software]][row,][[paste0(software,"_Individual_Assignment")]]), (results_list[[software]][[paste0(software,"_Individual_Assignment")]]))
					}
				} else{
					### Join to the datatable
					results_list[[software]] <- cluster_df[[software]][,.SD, .SD = c(paste0(software,"_Cluster"), paste0(software, "_Individual_Assignment"))][results_list[[software]], on = paste0(software,"_Cluster")]
				}
			} else {
				message("None of the ", software, " clusters had at least ", args$pct_agreement, " with an individual/cluster in ", ref)
			}

			### Update doublets 
			results_list[[software]][[paste0(software, "_Individual_Assignment")]] <- ifelse(is.na(results_list[[software]][[paste0(software, "_Individual_Assignment")]]), "doublet", results_list[[software]][[paste0(software, "_Individual_Assignment")]])
			results_list[[software]][[paste0(software, "_DropletType")]] <- ifelse(results_list[[software]][[paste0(software, "_Individual_Assignment")]] == "doublet", "doublet", results_list[[software]][[paste0(software, "_DropletType")]])
		}
		combined_results <- reduce(results_list, left_join)

	} else if (is.null(args$demuxlet) & is.null(args$vireo) & length(cluster_no_assign) == 1){
		### If only had one software => make clusters into assignments if assignments aren't there
		results_list[[cluster_no_assign[1]]][[paste0(cluster_no_assign[1], "_Individual_Assignment")]] <-  results_list[[cluster_no_assign[1]]][[paste0(cluster_no_assign[1], "_Cluster")]]
		combined_results <- reduce(results_list, left_join)
	
	}


	### Update any missing data with doublet ###
	combined_results[is.na(combined_results)] <- "unassigned"

	message("\nWriting output.\n")
	fwrite(combined_results, args$out, sep = "\t", append = FALSE)

	message("\nMaking and writing summary files.\n")
	combined_results_summary <- combined_results[,.(N=.N), by = c(grep("DropletType",colnames(combined_results), value = TRUE))]

	if (!is.null(args$demuxlet) | !is.null(args$freemuxlet) | !is.null(args$scSplit) | !is.null(args$souporcell) | !is.null(args$vireo)){
		columns <- grep("Individual_Assignment", colnames(combined_results), value = TRUE)
		if (!("Freemuxlet_Individual_Assignment" %in% columns) & "Freemuxlet_Cluster" %in% colnames(combined_results)){
			columns <- c(columns, "Freemuxlet_Cluster")
		}
		if (!("scSplit_Individual_Assignment" %in% columns) & "scSplit_Cluster" %in% colnames(combined_results)){
			columns <- c(columns, "scSplit_Cluster")
		}
		if (!("souporcell_Individual_Assignment" %in% columns) & "souporcell_Cluster" %in% colnames(combined_results)){
			columns <- c(columns, "souporcell_Cluster")
		}
		print(columns)
		demultiplex_combined_results_summary <- combined_results[,.(N=.N), by = c(columns)]
		demultiplex_combined_results_summary <- demultiplex_combined_results_summary[order(-N)]
		fwrite(demultiplex_combined_results_summary, paste0(tools::file_path_sans_ext(args$out),"_demultiplexing_summary.tsv"), sep = "\t", append = FALSE)

	}
	combined_results_summary <- combined_results_summary[order(-N)]

	fwrite(combined_results_summary, paste0(tools::file_path_sans_ext(args$out),"_summary.tsv"), sep = "\t", append = FALSE)
	

	## Call singlets and doublets using combined intersectional methods ##
	if (is.null(args$method)){
		message("No intersection classification method was provided. Output does not contain combined software cell type or donor classifications.")
	} else {
		if (!(args$method %in% c("MajoritySinglet", "AtLeastHalfSinglet", "AnySinglet", "AnyDoublet"))){
			message("Did not recognize the specified method. Can be one of 'MajoritySinglet', 'AtLeastHalfSinglet', 'AnySinglet' or 'AnyDoublet'. Output does not contain combined software cell type or donor classifications.")
		} else{
			### Check if using demultiplexing softwares and pull the individual ids if doing so
			if(any(grepl("Individual_Assignment", colnames(combined_results)))){

				individual_assignment_list <- future_apply(combined_results[,.SD, .SDcols = grep("Individual_Assignment", colnames(combined_results), value = TRUE)], 1, function(y) table(y)[!(rownames(table(y)) %in% c("unassigned","doublet"))])

				individual_assignment <- data.table(ID = unlist(lapply(individual_assignment_list, function(y) ifelse(is.null(names(which.max.simple(y[names(y) != "doublet"]))), NA, names(which.max.simple(y[names(y) != "doublet"]))))),
													N =  unlist(lapply(individual_assignment_list, function(y) ifelse(is.null(names(which.max.simple(y[names(y) != "doublet"]))), NA, max(y[names(y) != "doublet"]))))) ## Need to remove doublet count from table
			}

			### Make combined calls dependent on the method selected by the user
			if (args$method == "MajoritySinglet"){ ### Only call a singlet if majority of softwares call a singlet AND majority of demultiplexing (if present) call the same donor
				message("Using the 'MajoritySinglet' method to call droplet classification (and donor identity if demultiplexing softwares included).")
				
				if (!is.null(individual_assignment)){
					### method when have demultiplexing softwares
					combined_results$MajoritySinglet_DropletType <- ifelse(rowSums(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)] == "singlet") > length(grep("DropletType", colnames(combined_results)))/2 , "singlet", "doublet")
					combined_results$MajoritySinglet_Individual_Assignment <- ifelse(combined_results$MajoritySinglet_DropletType == "singlet" & !is.na(individual_assignment$ID), individual_assignment$ID,"doublet")
					combined_results$MajoritySinglet_DropletType <- ifelse(is.na(combined_results$MajoritySinglet_DropletType) & combined_results$MajoritySinglet_Individual_Assignment == "doublet", "doublet", combined_results$MajoritySinglet_DropletType)
					combined_results$MajoritySinglet_Individual_Assignment <- ifelse(combined_results$MajoritySinglet_DropletType == "singlet" & combined_results$MajoritySinglet_Individual_Assignment == "doublet", "unassigned", combined_results$MajoritySinglet_Individual_Assignment)
				} else {
					### method when no demultiplexing softwares
					combined_results$MajoritySinglet_DropletType <- ifelse(rowSums(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)] == "singlet") > length(grep("DropletType", colnames(combined_results)))/2, "singlet", "doublet")
				}
			} else if (args$method == "AtLeastHalfSinglet"){ ### Only call a singlet if at least half of softwares call a singlet AND at least half of demultiplexing (if present) call the same donor
				message("Using the 'AtLeastHalfSinglet' method to call droplet classification (and donor identify if demultiplexing softwares included).")
			
				if (any(!is.null(individual_assignment))){
					### method when have demultiplexing softwares
					combined_results$AtLeastHalfSinglet_DropletType <- ifelse(rowSums(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)] == "singlet") >= length(grep("DropletType", colnames(combined_results)))/2 , "singlet", "doublet")
					combined_results$AtLeastHalfSinglet_Individual_Assignment <- ifelse(combined_results$AtLeastHalfSinglet_DropletType == "singlet" & !is.na(individual_assignment$ID), individual_assignment$ID, "doublet")
					combined_results$AtLeastHalfSinglet_DropletType <- ifelse((is.na(combined_results$AtLeastHalfSinglet_DropletType) & combined_results$AtLeastHalfSinglet_Individual_Assignment == "doublet"), "doublet", combined_results$AtLeastHalfSinglet_DropletType)
					combined_results$AtLeastHalfSinglet_Individual_Assignment <- ifelse(combined_results$AtLeastHalfSinglet_DropletType == "singlet" & combined_results$AtLeastHalfSinglet_Individual_Assignment == "doublet", "unassigned", combined_results$AtLeastHalfSinglet_Individual_Assignment)
				} else {
					### method when no demultiplexing softwares
					combined_results$AtLeastHalfSinglet_DropletType <- ifelse(rowSums(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)] == "singlet") >= length(grep("DropletType", colnames(combined_results)))/2, "singlet", "doublet")
				}
			} else if (args$method == "AnySinglet"){ ### Call a singlet if any softwares calls that droplet a singlet AND using the donor with the most consensus (call doublet if no consensus)
				message("Using the 'AnySinglet' method to call droplet classification (and donor identify if demultiplexing softwares included).")
			
				if (any(!is.null(individual_assignment))){
					### method when have demultiplexing softwares
					combined_results$AnySinglet_DropletType <- ifelse(rowSums(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)] == "singlet") > 0 & !is.na(individual_assignment$N), "singlet", "doublet")
					combined_results$AnySinglet_Individual_Assignment <- ifelse(combined_results$AnySinglet_DropletType == "singlet" & !is.na(individual_assignment$ID), individual_assignment$ID, "doublet")
					combined_results$AnySinglet_DropletType <- ifelse(is.na(combined_results$AnySinglet_DropletType) & combined_results$AnySinglet_Individual_Assignment == "doublet", "doublet", combined_results$AnySinglet_DropletType)
					combined_results$AnySinglet_Individual_Assignment <- ifelse(combined_results$AnySinglet_DropletType == "singlet" & combined_results$AnySinglet_Individual_Assignment == "doublet", "unassigned", combined_results$AnySinglet_Individual_Assignment)
				} else {
				### method when no demultiplexing softwares
					combined_results$AnySinglet_DropletType <- ifelse(rowSums(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)] == "singlet") > 0, "singlet", "doublet")
				}
			} else if (args$method == "AnyDoublet"){ ### Call a singlet if all softwares calls that droplet a singlet AND all softwares call the same donor
				message("Using the 'AnyDoublet' method to call droplet classification (and donor identify if demultiplexing softwares included).")
			
				if (any(!is.null(individual_assignment))){
					### method when have demultiplexing softwares
					combined_results$AnyDoublet_DropletType <- ifelse(rowSums(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)] != "singlet") > 0 | is.na(individual_assignment$N), "doublet", "singlet")
					combined_results$AnyDoublet_Individual_Assignment <- ifelse(combined_results$AnyDoublet_DropletType == "doublet" | is.na(individual_assignment$ID), "doublet", individual_assignment$ID)
					combined_results$AnyDoublet_DropletType <- ifelse(is.na(combined_results$AnyDoublet_DropletType) & combined_results$AnyDoublet_Individual_Assignment == "doublet", "doublet", combined_results$AnyDoublet_DropletType)
					combined_results$AnyDoublet_Individual_Assignment <- ifelse(combined_results$AnyDoublet_DropletType == "singlet" & combined_results$AnyDoublet_Individual_Assignment == "doublet", "unassigned", combined_results$AnyDoublet_Individual_Assignment)
				} else {
					### method when no demultiplexing softwares
					combined_results$AnyDoublet_DropletType <- ifelse(rowSums(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)] == "doublet") > 0, "doublet", "singlet")
				}
			}
		}
	message("\nWriting output with combined calls.\n")
	fwrite(combined_results, paste0(tools::file_path_sans_ext(args$out),"_w_combined_assignments.tsv"), sep = "\t", append = FALSE)
	}

	##### Make an upset plot for the results to visualize the agreement between different softwares #####
	### First check whether common assignments present (if they are, then will color by common assignment in the upset bar plots) ###
	if (any(colnames(combined_results) %in% c("AnyDoublet_Individual_Assignment", "AnySinglet_Individual_Assignment", "AtLeastHalfSinglet_Individual_Assignment", "MajoritySinglet_Individual_Assignment"))){
		## Dataframe with each software singlet (1), doublet (0) designations + column for individual assignment
		upset_df <- data.frame(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)])
		upset_df <- upset_df[, !(colnames(upset_df) %in% c("MajoritySinglet_DropletType", "AtLeastHalfSinglet", "AnySinglet_DropletType", "AnyDoublet_DropletType"))]
		upset_df <- as.data.frame(1*(upset_df == "singlet"))
		colnames(upset_df) <- gsub("_DropletType", "", colnames(upset_df))
		softwares_ordered <- rev(c("Demuxlet", "Freemuxlet", "scSplit", "Souporcell", "Vireo", "DoubletDecon", "DoubletDetection", "DoubletFinder", "scDblFinder", "scds", "scrublet", "solo"))
		columns <- softwares_ordered[softwares_ordered %in% colnames(upset_df)]
		upset_df <- upset_df[,c(columns)]

		upset_df$Final_Individual_Assignment <- as.vector(combined_results[,.SD, .SDcols = (colnames(combined_results) %in% c("AnySinglet_Individual_Assignment", "AtLeastHalfSinglet_Individual_Assignment", "AnyDoublet_Individual_Assignment", "MajoritySinglet_Individual_Assignment"))][[1]])
		upset_df$Final_Individual_Assignment <- factor(upset_df$Final_Individual_Assignment, c(sort(unique(upset_df$Final_Individual_Assignment)[!(unique(upset_df$Final_Individual_Assignment) %in% c("doublet", "unassigned"))]), "doublet", "unassigned"))


		colourCount = length(unique(upset_df$Final_Individual_Assignment)[!unique(upset_df$Final_Individual_Assignment) %in% c("unassigned", "doublet")])
		getPalette = colorRampPalette(c("#f44336", "#e81e63", "#9c27b0", "#673ab7", "#3f51b5", "#2196f3", "#03a9f4", "#00bcd4", "#009688", "#4caf50", "#8bc34a", "#cddc39", "#ffeb3b","#ffc107", "#ff9800", "#ff5722", "#795548"))

		colors <- getPalette(colourCount)
		names(colors) <- unique(upset_df$Final_Individual_Assignment)[!unique(upset_df$Final_Individual_Assignment) %in% c("unassigned", "doublet")]

		if ("doublet" %in% unique(upset_df$Final_Individual_Assignment)){
			colors <- c(colors, doublet = "#000000")
		}
		if ("unassigned" %in% unique(upset_df$Final_Individual_Assignment)){
			colors <- c(colors, unassigned = "grey70")
		}

		pUpset <- upset(upset_df,
			columns,
			# set_sizes=TRUE, 
			sort_sets=FALSE,
			name="Singlet Classifications",
			set_sizes=(
				upset_set_size()
				+ ylab('Number Singlets')
				+ theme(axis.text.x=element_text(angle=90))
    		),
			base_annotations=list(
				'Intersection size'=intersection_size(
					counts=FALSE,
					mapping=aes(fill=Final_Individual_Assignment)
				) +
					ylab('Number Droplets') +
					scale_fill_manual(values = colors, name = "Final\nIndividual\nAssignment")
			),
			width_ratio=0.1
		)

		message("\nMaking figures of software and final assignments.\n")

		pdf(file= paste0(tools::file_path_sans_ext(args$out),"Singlets_upset.pdf"), height = 5, width = 10) # or other device
		print(pUpset)
		dev.off()


	} else if (any(colnames(combined_results) %in% c("MajoritySinglet_DropletType", "AtLeastHalfSinglet_DropletType", "AnySinglet_DropletType", "AnyDoublet_DropletType"))){
		## Dataframe with each software singlet (1), doublet (0) designations
		upset_df <- data.frame(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)])
		upset_df <- upset_df[, !(colnames(upset_df) %in% c("MajoritySinglet_DropletType", "AtLeastHalfSinglet_DropletType", "AnySinglet_DropletType", "AnyDoublet_DropletType"))]
		upset_df <- as.data.frame(1*(upset_df == "singlet"))
		colnames(upset_df) <- gsub("_DropletType", "", colnames(upset_df))
		softwares_ordered <- rev(c("Demuxlet", "Freemuxlet", "scSplit", "Souporcell", "Vireo", "DoubletDecon", "DoubletDetection", "DoubletFinder", "scDblFinder", "scds", "scrublet", "solo"))
		columns <- softwares_ordered[softwares_ordered %in% colnames(upset_df)]
		upset_df <- upset_df[,c(columns)]

		upset_df$Final_Assignment <- as.vector(combined_results[,.SD, .SDcols = (colnames(combined_results) %in% c("MajoritySinglet_DropletType", "AtLeastHalfSinglet_DropletType", "AnySinglet_DropletType", "AnyDoublet_DropletType"))][[1]])
		upset_df$Final_Assignment <- factor(upset_df$Final_Assignment, c(sort(unique(upset_df$Final_Assignment)[!(unique(upset_df$Final_Assignment) %in% c("doublet", "unassigned"))]), "doublet", "unassigned"))


		colourCount = length(unique(upset_df$Final_Assignment))
		getPalette = colorRampPalette(c("#f44336", "#e81e63", "#9c27b0", "#673ab7", "#3f51b5", "#2196f3", "#03a9f4", "#00bcd4", "#009688", "#4caf50", "#8bc34a", "#cddc39", "#ffeb3b","#ffc107", "#ff9800", "#ff5722", "#795548", "#9e9e9e", "#607d8b", "#000000",))

		pUpset <- upset(upset_df,
			columns,
			sort_sets=FALSE,
			name="Singlet Classifications",
			set_sizes=(
				upset_set_size()
				+ ylab('Number Singlets')
				+ theme(axis.text.x=element_text(angle=90))
    		),
			base_annotations=list(
				'Intersection size'=intersection_size(
					counts=FALSE,
					mapping=aes(fill=Final_Assignment)
				) +
					ylab('Number Droplets') +
					scale_fill_manual(values = getPalette(colourCount), name = "Final\nAssignment")
			),
			width_ratio=0.1
		)
		
		message("\nMaking figures of software and final assignments.\n")

		pdf(file= paste0(tools::file_path_sans_ext(args$out),"Singlets_upset.pdf"), height = 5, width = 10) # or other device
		print(pUpset)
		dev.off()

	} else {
		## Dataframe with each software singlet (1), doublet (0) designations
		upset_df <- data.frame(combined_results[,.SD, .SDcols = grep("DropletType", colnames(combined_results), value = TRUE)])
		upset_df <- upset_df[, !(colnames(upset_df) %in% c("MajoritySinglet_DropletType", "AtLeastHalfSinglet_DropletType", "AnySinglet_DropletType", "AnyDoublet_DropletType"))]
		upset_df <- as.data.frame(1*(upset_df == "singlet"))
		colnames(upset_df) <- gsub("_DropletType", "", colnames(upset_df))
		softwares_ordered <- rev(c("Demuxlet", "Freemuxlet", "scSplit", "Souporcell", "Vireo", "DoubletDecon", "DoubletDetection", "DoubletFinder", "scDblFinder", "scds", "scrublet", "solo"))
		columns <- softwares_ordered[softwares_ordered %in% colnames(upset_df)]
		upset_df <- upset_df[,c(columns)]

		pUpset <- upset(upset_df,
			columns,
			sort_sets=FALSE,
			name="Singlet Classifications",
			set_sizes=(
				upset_set_size()
				+ ylab('Number Singlets')
				+ theme(axis.text.x=element_text(angle=90))
    		),
			base_annotations=list(
				'Intersection size'=intersection_size(
					counts=FALSE
				) +
					ylab('Number Droplets') 			),
			width_ratio=0.1
		)

		message("\nMaking figures of final assignments.\n")

		pdf(file= paste0(tools::file_path_sans_ext(args$out),"Singlets_upset.pdf"), height = 5, width = 10) # or other device
		print(pUpset)
		dev.off()

	}


	message("\nDone!\n")
}
