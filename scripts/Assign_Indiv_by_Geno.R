#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")
library(argparse)


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-r", "--reference_vcf", required = TRUE, type = "character", help="The output directory where results will be saved")
parser$add_argument("-c", "--cluster_vcf", required = TRUE, type = "character", help = "A QC, normalized seurat object with classifications/clusters as Idents().")
parser$add_argument("-o", "--outdir", required = TRUE, type = "character", help = "Number  of genes to use in \'Improved_Seurat_Pre_Process\' function.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()


library(tidyr)
library(tidyverse)
library(dplyr)
library(vcfR)
library(lsa)
library(ComplexHeatmap)

########## Set up functions ##########
calculate_DS <- function(GP_df){
    columns <- c()
    for (i in 1:ncol(GP_df)){
        columns <- c(columns, paste0(colnames(GP_df)[i],"-0"), paste0(colnames(GP_df)[i],"-1"), paste0(colnames(GP_df)[i],"-2"))
    }
    df <- GP_df
    colnames(df) <- paste0("c", colnames(df))
    colnames_orig <- colnames(df)
    for (i in 1:length(colnames_orig)){
        df <- separate(df, sep = ",", col = colnames_orig[i], into = columns[(1+(3*(i-1))):(3+(3*(i-1)))])
    }
    df <- mutate_all(df, function(x) as.numeric(as.character(x)))
    for (i in 1: ncol(GP_df)){
        GP_df[,i] <- df[,(2+((i-1)*3))] + 2* df[,(3+((i-1)*3))]
    }
    return(GP_df)
}

pearson_correlation <- function(df, ref_df, clust_df){
    for (col in colnames(df)){
        for (row in rownames(df)){
            df[row,col] <- cor(as.numeric(pull(ref_df, col)), as.numeric(pull(clust_df, row)), method = "pearson", use = "complete.obs")
        }
    }
    return(df)
}


########## Read in vcf files for each of three non-reference genotype softwares ##########
ref_geno <- read.vcfR(args$reference_vcf)
cluster_geno <- read.vcfR(args$cluster_vcf)



########## Convert to tidy data frame ##########
####### Identify which genotype FORMAT to use #######
##### Cluster VCF #####
### Check for each of the different genotype formats ##
## DS ##
format_clust=NA
cluster_geno_tidy <- as_tibble(extract.gt(element = "DS",cluster_geno, IDtoRowNames = F))
if (!all(colSums(is.na(cluster_geno_tidy)) == nrow(cluster_geno_tidy))){
	message("Found DS genotype format in cluster vcf. Will use that metric for cluster correlation.")
	format_clust = "DS"
}

## GT ##
if (is.na(format_clust)){
	cluster_geno_tidy <- as_tibble(extract.gt(element = "GT",cluster_geno, IDtoRowNames = F))
	if (!all(colSums(is.na(cluster_geno_tidy)) == nrow(cluster_geno_tidy))){
		message("Found GT genotype format in cluster vcf. Will use that metric for cluster correlation.")
		format_clust = "GT"

		if (any(grepl("\\|",cluster_geno_tidy[1,]))){
			separator = "|"
			message("Detected | separator for GT genotype format in cluster vcf")
		} else if (any(grepl("/",cluster_geno_tidy[1,]))) {
			separator = "/"
			message("Detected / separator for GT genotype format in cluster vcf")
		} else {
			format_clust = NA
			message("Can't identify a separator for the GT field in cluster vcf, moving on to using GP.")
		}

		cluster_geno_tidy <- as_tibble(lapply(cluster_geno_tidy, function(x) {gsub(paste0("0",separator,"0"),0, x)}) %>%
		                        lapply(., function(x) {gsub(paste0("0",separator,"1"),1, x)}) %>%
		                        lapply(., function(x) {gsub(paste0("1",separator,"0"),1, x)}) %>%
		                        lapply(., function(x) {gsub(paste0("1",separator,"1"),2, x)}))

	}
}

## GP ##
if (is.na(format_clust)){
	cluster_geno_tidy <- as_tibble(extract.gt(element = "GP",cluster_geno, IDtoRowNames =F))
	if (!all(colSums(is.na(cluster_geno_tidy)) == nrow(cluster_geno_tidy))){
		format_clust = "GP"
		cluster_geno_tidy <- calculate_DS(cluster_geno_tidy)
		message("Found GP genotype format in cluster vcf. Will use that metric for cluster correlation.")

	} else {
		print("Could not identify the expected genotype format fields (DS, GT or GP) in your cluster vcf. Please check the vcf file and make sure that one of the expected genotype format fields is included or run manually with your genotype format field of choice. Quitting")
		q()
	}
}

    



### Reference VCF ###
### Check for each of the different genotype formats ##
## DS ##
format_ref = NA
ref_geno_tidy <- as_tibble(extract.gt(element = "DS",ref_geno, IDtoRowNames = F))
if (!all(colSums(is.na(ref_geno_tidy)) == nrow(ref_geno_tidy))){
	message("Found DS genotype format in reference vcf. Will use that metric for cluster correlation.")
	format_ref = "DS"
}

## GT ##
if (is.na(format_ref)){
	ref_geno_tidy <- as_tibble(extract.gt(element = "GT",ref_geno, IDtoRowNames = F))
	if (!all(colSums(is.na(ref_geno_tidy)) == nrow(ref_geno_tidy))){
		message("Found GT genotype format in reference vcf. Will use that metric for cluster correlation.")
		format_ref = "GT"

		if (any(grepl("\\|",ref_geno_tidy[1,]))){
			separator = "|"
			message("Detected | separator for GT genotype format in reference vcf")
		} else if (any(grepl("/",ref_geno_tidy[1,]))) {
			separator = "/"
			message("Detected / separator for GT genotype format in reference vcf")
		} else {
			format_ref = NA
			message("Can't identify a separator for the GT field in reference vcf, moving on to using GP.")
		}

		ref_geno_tidy <- as_tibble(lapply(ref_geno_tidy, function(x) {gsub(paste0("0",separator,"0"),0, x)}) %>%
		                        lapply(., function(x) {gsub(paste0("0",separator,"1"),1, x)}) %>%
		                        lapply(., function(x) {gsub(paste0("1",separator,"0"),1, x)}) %>%
		                        lapply(., function(x) {gsub(paste0("1",separator,"1"),2, x)}))

	}
}

## GP ##
if (is.na(format_ref)){
	ref_geno_tidy <- as_tibble(extract.gt(element = "GP",ref_geno, IDtoRowNames = F))
	if (!all(colSums(is.na(ref_geno_tidy)) == nrow(ref_geno_tidy))){
		format_clust = "GP"
		ref_geno_tidy <- calculate_DS(ref_geno_tidy)
		message("Found GP genotype format in cluster vcf. Will use that metric for cluster correlation.")

	} else {
		print("Could not identify the expected genotype format fields (DS, GT or GP) in your cluster vcf. Please check the vcf file and make sure that one of the expected genotype format fields is included or run manually with your genotype format field of choice. Quitting")
		q()
	}
}



### Get SNP IDs that will match between reference and cluster ###
## Account for possibility that the ref or alt might be missing
if ((all(is.na(cluster_geno@fix[,'REF'])) & all(is.na(cluster_geno@fix[,'ALT']))) | (all(is.na(ref_geno@fix[,'REF'])) & all(is.na(ref_geno@fix[,'ALT'])))){
	message("The REF and ALT categories are not provided for the reference and/or the cluster vcf. Will use just the chromosome and position to match SNPs.")
	cluster_geno_tidy$ID <- paste0(cluster_geno@fix[,'CHROM'],":", cluster_geno@fix[,'POS'])
	ref_geno_tidy$ID <- paste0(ref_geno@fix[,'CHROM'],":", ref_geno@fix[,'POS'])
} else if (all(is.na(cluster_geno@fix[,'REF'])) | all(is.na(ref_geno@fix[,'REF']))){
	message("The REF categories are not provided for the reference and/or the cluster vcf. Will use the chromosome, position and ALT to match SNPs.")
	cluster_geno_tidy$ID <- paste0(cluster_geno@fix[,'CHROM'],":", cluster_geno@fix[,'POS'],"_", cluster_geno@fix[,'REF'])
	ref_geno_tidy$ID <- paste0(ref_geno@fix[,'CHROM'],":", ref_geno@fix[,'POS'],"_", ref_geno@fix[,'REF'])
} else if (all(is.na(cluster_geno@fix[,'ALT'])) | all(is.na(ref_geno@fix[,'ALT']))){
	message("The ALT categories are not provided for the reference and/or the cluster vcf. Will use the chromosome, position and REF to match SNPs.")
	cluster_geno_tidy$ID <- paste0(cluster_geno@fix[,'CHROM'],":", cluster_geno@fix[,'POS'],"_", cluster_geno@fix[,'ALT'])
	ref_geno_tidy$ID <- paste0(ref_geno@fix[,'CHROM'],":", ref_geno@fix[,'POS'],"_", ref_geno@fix[,'ALT'])
} else {
	message("Found REF and ALT in both cluster and reference genotype vcfs. Will use chromosome, position, REF and ALT to match SNPs.")
		cluster_geno_tidy$ID <- paste0(cluster_geno@fix[,'CHROM'],":", cluster_geno@fix[,'POS'],"_", cluster_geno@fix[,'REF'],"_", cluster_geno@fix[,'ALT'])
	ref_geno_tidy$ID <- paste0(ref_geno@fix[,'CHROM'],":", ref_geno@fix[,'POS'],"_", ref_geno@fix[,'REF'],"_", ref_geno@fix[,'ALT'])
}


### Update the vcf dfs to remove SNPs with no genotyopes
cluster_geno_tidy <- cluster_geno_tidy[colSums(!is.na(cluster_geno_tidy)) > 0]
ref_geno_tidy <- ref_geno_tidy[colSums(!is.na(ref_geno_tidy)) > 0]



########## Get a unique list of SNPs that is in both the reference and cluster genotypes ##########
locations  <- inner_join(ref_geno_tidy[,"ID"],cluster_geno_tidy[,"ID"])
locations <- locations[!(locations$ID %in% locations[duplicated(locations),]$ID),]

########## Keep just the SNPs that overlap ##########
ref_geno_tidy <- left_join(locations, ref_geno_tidy)
cluster_geno_tidy <- left_join(locations, cluster_geno_tidy)

########## Correlate all the cluster genotypes with the individuals genotyped ##########
##### Make a dataframe that has the clusters as the row names and the individuals as the column names #####
pearson_correlations <- as.data.frame(matrix(nrow = (ncol(cluster_geno_tidy) -1), ncol = (ncol(ref_geno_tidy) -1)))
colnames(pearson_correlations) <- colnames(ref_geno_tidy)[2:(ncol(ref_geno_tidy))]
rownames(pearson_correlations) <- colnames(cluster_geno_tidy)[2:(ncol(cluster_geno_tidy))]
pearson_correlations <- pearson_correlation(pearson_correlations, ref_geno_tidy, cluster_geno_tidy)
cluster <- data.frame("Cluster" = rownames(pearson_correlations))
pearson_correlations_out <- cbind(cluster, pearson_correlations)

########## Save the correlation dataframes ##########
write_delim(pearson_correlations_out, file = paste0(args$outdir,"/ref_clust_pearson_correlations.tsv"), delim = "\t" )


########## Create correlation figures ##########
col_fun = colorRampPalette(c("white", "red"))(101)
pPearsonCorrelations <- Heatmap(as.matrix(pearson_correlations), cluster_rows = T, col = col_fun)

########## Save the correlation figures ##########
png(filename = paste0(args$outdir,"/ref_clust_pearson_correlation.png"), width = 500)
print(pPearsonCorrelations)
dev.off()

########## Assign individual to cluster based on highest correlating individual ##########
key <- as.data.frame(matrix(nrow = ncol(pearson_correlations), ncol = 3))
colnames(key) <- c("Genotype_ID","Cluster_ID","Correlation")
key$Genotype_ID <- colnames(pearson_correlations)
for (id in key$Genotype_ID){
    if (max(pearson_correlations[,id]) == max(pearson_correlations[rownames(pearson_correlations)[which.max(pearson_correlations[,id])],])){
        key$Cluster_ID[which(key$Genotype_ID == id)] <- rownames(pearson_correlations)[which.max(pearson_correlations[,id])]
        key$Correlation[which(key$Genotype_ID == id)] <- max(pearson_correlations[,id])
    } else {
        key$Cluster_ID[which(key$Genotype_ID == id)] <- "unassigned"
        key$Correlation[which(key$Genotype_ID == id)] <- NA
    }
}

write_delim(key, file =paste0(args$outdir,"/Genotype_ID_key.txt"), delim = "\t")


