#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(tidyverse)))



# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, type = "character", help="The file where results will be saved")
parser$add_argument("-d", "--demuxlet", required = FALSE, type = "character", default = NULL, help = "Path to demuxlet results. Only use this option if you want to include the demuxlet results.")
parser$add_argument("-f", "--freemuxlet", required = FALSE, type = "character", default=NULL, help = "Path to freemuxlet results. Only use this option if you want to include the freemuxlet results.")
parser$add_argument("-g", "--freemuxlet_assignments", required = FALSE, type = "character", default=NULL, help = "Path to freemuxlet cluster-to-individual assignments. Only use this option if have used reference SNP genotypes to assign individuals to clusters for the freemuxlet results.")
parser$add_argument("-s", "--scSplit", required = FALSE, type="character", default=NULL, help="Path to scSplit results. Only use this option if you want to include the scSplit results.")
parser$add_argument("-w", "--scSplit_assignments", required = FALSE, type="character", default=NULL, help="Path to scSplit cluster-to-individual assignments. Only use this option if you have used reference SNP genotypes to assign individuals to clusters for the scSplit results.")
parser$add_argument("-u", "--souporcell", required = FALSE, type = "character", default=NULL, help = "Path to souporcell results. Only use this option if you want to include the souporcell results.")
parser$add_argument("-x", "--souporcell_assignments", required = FALSE, type = "character", default=NULL, help = "Path to souporcell cluster-to-individual assignments. Only use this option if you have used reference SNP genotypes to assign individuals to clusters for the souporcell results.")
parser$add_argument("-v", "--vireo", required = FALSE, type = "character", default=NULL, help = "Path to vireo results. Only use this option if you want to include the vireo results.")
parser$add_argument("-e", "--DoubletDecon", required = FALSE, type = "character", default=NULL, help = "Path to DoubletDecon results. Only use this option if you want to include the DoubletDecon results.")
parser$add_argument("-t", "--DoubletDetection", required = FALSE, type = "character", default=NULL, help = "Path to DoubletDetection results. Only use this option if you want to include the DoubletDetection results.")
parser$add_argument("-i", "--DoubletFinder", required = FALSE, type = "character", default=NULL, help = "Path to DoubletFinder results. Only use this option if you want to include the DoubletFinder results.")
parser$add_argument("-n", "--scDblFinder", required = FALSE, type = "character", default=NULL, help = "Path to scDblFinder results. Only use this option if you want to include the scDblFinder results.")
parser$add_argument("-c", "--scds", required = FALSE, type = "character", default=NULL, help = "Path to scds results. Only use this option if you want to include the scds results.")
parser$add_argument("-r", "--scrublet", required = FALSE, type = "character", default=NULL, help = "Path to scrublet results. Only use this option if you want to include the scrublet results.")
parser$add_argument("-l", "--solo", required = FALSE, type = "character", default=NULL, help = "Path to solo results. Only use this option if you want to include the solo results.")
# parser$add_argument("-m", "--method", required = FALSE, type = "character", default=NULL, help = "Combination method. Leave blank if you just want all the combine the methods into a single dataframe.")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()



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
		results_list[["demuxlet"]]$Demuxlet_DropletType  <- gsub("SNG", "singlet", results_list[["demuxlet"]]$Demuxlet_DropletType) %>% gsub("DBL", "doublet", .) %>% gsub("AMB", "ambiguous", .)

		### Update Individual column to be individual or doublet or ambiguous/unassigned ###
		results_list[["demuxlet"]]$Demuxlet_Individual_Assignment <- ifelse(results_list[["demuxlet"]]$Demuxlet_DropletType  == "doublet", "doublet",
																		ifelse(results_list[["demuxlet"]]$Demuxlet_DropletType  == "unassigned", "unassigned", gsub(",.+,.\\..+","",results_list[["demuxlet"]]$Demuxlet_Individual_Assignment)))

	}


	if (!is.null(args$freemuxlet)){
		message("Reading in freemuxlet results.")
		if (file_test("-f", args$freemuxlet)){
			results_list[["freemuxlet"]] <- fread(args$freemuxlet, sep = "\t")
		} else {
			freemuxlet <- list.files(args$freemuxlet, pattern = ".samples.gz")
			if (length(freemuxlet) == 1){
					results_list[["freemuxlet"]] <- fread(paste0(args$freemuxlet, "/", freemuxlet), sep = "\t")
			} else {
				message("Can't read in the freemuxlet report either from the file you provided or to find the '*.samples.gz' file in the directory you provided. Please double check your path and provide the full path to the freemuxlet file. Exiting.")
				q()
			}
		}

		### Update dataframe to just be categories of interest
		results_list[["freemuxlet"]] <- results_list[["freemuxlet"]][,c("BARCODE", "DROPLET.TYPE", "BEST.GUESS")]
		colnames(results_list[["freemuxlet"]]) <- c("Barcode", "Freemuxlet_DropletType", "Freemuxlet_Cluster")

		### update Demuxlet_DropletType  column SNG and DBL ###
		results_list[["freemuxlet"]]$Freemuxlet_DropletType  <- gsub("SNG", "singlet", results_list[["freemuxlet"]]$Freemuxlet_DropletType ) %>% gsub("DBL", "doublet", .) %>% gsub("AMB", "ambiguous", .)

		### Update Individual column to be individual or doublet or ambiguous/unassigned ###
		results_list[["freemuxlet"]]$Freemuxlet_Cluster <- ifelse(results_list[["freemuxlet"]]$Freemuxlet_DropletType  == "doublet", "doublet",
																		ifelse(results_list[["freemuxlet"]]$Freemuxlet_DropletType  == "unassigned", "unassigned", gsub(",.+","",results_list[["freemuxlet"]]$Freemuxlet_Cluster)))

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
			results_list[["souporcell"]] <- fread(args$souporcell, sep = "\t")
		} else {
			tryCatch({
				results_list[["souporcell"]] <- fread(paste0(args$souporcell, "/clusters.tsv"), sep = "\t")
			},
			error = function(e) {
				message(e)
				message("Can't read in the souporcell report either from the file you provided or to find the 'clusters.tsv' file in the directory you provided. Please double check your path and provide the full path to the souporcell file. Exiting.")
				q()
			})
		}

		### Update dataframe to just be categories of interest
		results_list[["souporcell"]] <- results_list[["souporcell"]][,c("barcode", "status","assignment")]
		colnames(results_list[["souporcell"]]) <- c("Barcode", "Souporcell_DropletType", "Souporcell_Cluster")
		results_list[["souporcell"]]$Souporcell_Cluster <- ifelse(results_list[["souporcell"]]$Souporcell_DropletType  == "doublet", "doublet", results_list[["souporcell"]]$Souporcell_Cluster)

	}

	if (!is.null(args$vireo)){
		message("Reading in vireo results.")
		if (file_test("-f", args$vireo)){
			results_list[["vireo"]] <- fread(args$vireo, sep = "\t")
		} else {
			tryCatch({
				results_list[["vireo"]] <- fread(paste0(args$vireo, "/donor_ids.tsv"), sep = "\t")
			},
			error = function(e) {
				message(e)
				message("Can't read in the vireo report either from the file you provided or to find the 'donor_ids.tsv' file in the directory you provided. Please double check your path and provide the full path to the vireo file. Exiting.")
				q()
			})
		}
		
		### Update dataframe to just be categories of interest
		results_list[["vireo"]] <- results_list[["vireo"]][,c("cell", "donor_id")]
		colnames(results_list[["vireo"]]) <- c("Barcode", "vireo_Individual_Assignment")
		results_list[["vireo"]]$vireo_DropletType  <- ifelse(results_list[["vireo"]]$vireo_Individual_Assignment == "doublet", "doublet", "singlet")

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
				results_assignments_list[["freemuxlet"]] <- fread(args$freemuxlet_assignments, sep = "\t")
			} else {
				tryCatch({
					results_assignments_list[["freemuxlet"]] <- fread(paste0(args$freemuxlet_assignments, "/Genotype_ID_key.txt"), sep = "\t")
				},
				error = function(e) {
					message("Can't read in the freemuxlet cluster-to-indiviudal assignments report either from the file you provided or to find the 'Genotype_ID_key.tsv' file in the directory you provided. Please double check your path and provide the full path to the freemuxlet cluster-to-indiviudal key file. Exiting.")
					q()
				})
			}


			### Check for column names ###
			if (all(!(grepl("Genotype_ID",colnames(results_assignments_list[["freemuxlet"]]))) & !(grepl("Cluster_ID", colnames(results_assignments_list[["freemuxlet"]]))))){
				message("Didn't find 'Genotype_ID' and 'Cluster_ID' in the columns of your freemuxlet cluster-to-indiviudal assignments file. Will use first column as the individual ID and the second column will be the cluster ID.")
			}
			colnames(results_assignments_list[["freemuxlet"]])[1:2] <- c("Freemuxlet_Individual_Assignment", "Freemuxlet_Cluster")
			results_assignments_list[["freemuxlet"]]$Freemuxlet_Cluster <- as.character(results_assignments_list[["freemuxlet"]]$Freemuxlet_Cluster)

			### remove "CLUST" from assignments ###
			results_assignments_list[["freemuxlet"]]$Freemuxlet_Cluster <- gsub("CLUST", "", results_assignments_list[["freemuxlet"]]$Freemuxlet_Cluster)
			results_assignments_list[["freemuxlet"]]$Correlation <- NULL

			### Join
			message("Adding assignments to freeuxlet dataframe")
			results_list[["freemuxlet"]] <- results_assignments_list[["freemuxlet"]][results_list[["freemuxlet"]], on="Freemuxlet_Cluster"]  

			results_list[["freemuxlet"]]$Freemuxlet_Individual_Assignment <- ifelse(results_list[["freemuxlet"]]$Freemuxlet_DropletType == "doublet", "doublet", results_list[["freemuxlet"]]$Freemuxlet_Individual_Assignment)


			### check for excessive NA after joining ###
			if (length(which(is.na(results_list[["freemuxlet"]]$Freemuxlet_Individual_Assignment)))/nrow(results_list[["freemuxlet"]]) >= 0.5){
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
			}
			colnames(results_assignments_list[["scSplit"]])[1:2] <- c("scSplit_Individual_Assignment", "scSplit_Cluster")
			results_assignments_list[["scSplit"]]$scSplit_Cluster <- as.character(results_assignments_list[["scSplit"]]$scSplit_Cluster)

			results_assignments_list[["scSplit"]]$Correlation <- NULL


			### Join
			message("Adding assignments to scSplit dataframe")
			results_list[["scSplit"]] <- results_assignments_list[["scSplit"]][results_list[["scSplit"]], on="scSplit_Cluster"]  

			results_list[["scSplit"]]$scSplit_Individual_Assignment <- ifelse(results_list[["scSplit"]]$scSplit_DropletType == "doublet", "doublet", results_list[["scSplit"]]$scSplit_Individual_Assignment)

			### check for excessive NA after joining ###
			if (length(which(is.na(results_list[["scSplit"]]$Freemuxlet_Individual_Assignment)))/nrow(results_list[["scSplit"]]) >= 0.5){
				message("WARNING: There are more than 50% NA assignments after joining the cluster-to-individual assignments for scSplit. Please double check you files if this is unexpected.")
			}
		}

		if (!is.null(args$souporcell_assignments)){
			message("Reading in souporcell_assignments cluster-to-indiviudal assignments.")
			if (file_test("-f", args$souporcell_assignments)){
				results_assignments_list[["souporcell"]] <- fread(args$souporcell_assignments, sep = "\t")
			} else {
				tryCatch({
					results_assignments_list[["souporcell"]] <- fread(paste0(args$souporcell_assignments, "/Genotype_ID_key.txt"), sep = "\t")
				},
				error = function(e) {
					message("Can't read in the souporcell cluster-to-indiviudal assignments report either from the file you provided or to find the 'Genotype_ID_key.tsv' file in the directory you provided. Please double check your path and provide the full path to the souporcell cluster-to-indiviudal key file. Exiting.")
					q()
				})
			}


			### Check for column names ###
			if (all(!(grepl("Genotype_ID",colnames(results_assignments_list[["souporcell"]]))) & !(grepl("Cluster_ID", colnames(results_assignments_list[["souporcell"]]))))){
				message("Didn't find 'Genotype_ID' and 'Cluster_ID' in the columns of your souporcell cluster-to-indiviudal assignments. Will use first column as the individual ID and the second column will be the cluster ID.")
			}
			colnames(results_assignments_list[["souporcell"]])[1:2] <- c("Souporcell_Individual_Assignment", "Souporcell_Cluster")
			results_assignments_list[["souporcell"]]$Souporcell_Cluster <- as.character(results_assignments_list[["souporcell"]]$Souporcell_Cluster)
			results_assignments_list[["souporcell"]]$Correlation <- NULL


			### Join
			message("Adding assignments to souporcell dataframe")
			results_list[["souporcell"]] <- results_assignments_list[["souporcell"]][results_list[["souporcell"]], on="Souporcell_Cluster"]  

			results_list[["souporcell"]]$Souporcell_Individual_Assignment <- ifelse(results_list[["souporcell"]]$Souporcell_DropletType == "doublet", "doublet", results_list[["souporcell"]]$Souporcell_Individual_Assignment)

			### check for excessive NA after joining ###
			if (length(which(is.na(results_list[["souporcell"]]$Freemuxlet_Individual_Assignment)))/nrow(results_list[["souporcell"]]) >= 0.5){
				message("WARNING: There are more than 50% NA assignments after joining the cluster-to-individual assignments for souporcell. Please double check you files if this is unexpected.")
			}
			print(results_list[["souporcell"]] )
		}

	}

	message("\nCombining results.\n")
	combined_results <- reduce(results_list, left_join)

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
		fwrite(demultiplex_combined_results_summary, paste0(tools::file_path_sans_ext(args$out),"_demultiplexing_summary.tsv"), sep = "\t", append = FALSE)

	}
	fwrite(combined_results_summary, paste0(tools::file_path_sans_ext(args$out),"_summary.tsv"), sep = "\t", append = FALSE)
	

	message("\nDone!\n")
}


## IN docs: can provide file produced any way you want but must have either "Genotype_ID" and "Cluster_ID" as headers or the first must be individual ID and second column must be cluster id


## Possibly add methods to call singlets and doublets (probably in the future)


