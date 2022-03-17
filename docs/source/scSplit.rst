.. _scSplit-docs:

ScSplit
===========================

.. _ScSplit: https://github.com/jon-xu/scSplit
.. _preprint: https://www.biorxiv.org/content/10.1101/2022.03.07.483367v1

ScSplit_ is a reference-free demultiplexing software. If you have reference SNP genotypes, it would be better to use a demultiplexing software that can handle reference SNP genotypes (:ref:`Demuxlet <Demuxlet-docs>`, :ref:`Souporcell <Souporcell-docs>` or :ref:`Vireo<Vireo-docs>`)

Data
----
This is the data that you will need to have prepared to run ScSplit_:

.. admonition:: Required
  :class: important

  - Bam file (``$BAM``)

    - Aligned single cell reads

  - Genome reference fasta file (``$FASTA``)

  - Barcode file (``$BARCODES``)

  - Common SNP genotypes vcf (``$VCF``)

    - While not exactly required, using common SNP genotype locations enhances accuracy

      - If you have reference SNP genotypes for individuals in your pool, you can use those

      - If you do not have reference SNP genotypes, they can be from any large population resource (i.e. 1000 Genomes or HRC)

    - Filter for common SNPs (> 5% minor allele frequency) and SNPs overlapping genes

  - Number of samples in pool (``$N``)

  - Output directory (``$SCSPLIT_OUTDIR``)


Run ScSplit
-----------
Prepare Bam file
^^^^^^^^^^^^^^^^
First, you will need to prepare the bam file so that it only contains high quality, primarily mapped reads without any PCR duplicated reads.

.. code-block:: bash

  singularity exec Demuxafy.sif samtools view -b -S -q 10 -F 3844 $BAM > $SCSPLIT_OUTDIR/filtered_bam.bam
  singularity exec Demuxafy.sif samtools rmdup $SCSPLIT_OUTDIR/filtered_bam.bam $SCSPLIT_OUTDIR/filtered_bam_dedup.bam
  singularity exec Demuxafy.sif samtools sort -o $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam $SCSPLIT_OUTDIR/filtered_bam_dedup.bam
  singularity exec Demuxafy.sif samtools index $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam

After running these bam preparation steps, you will have the following files in your ``$SCSPLIT_OUTDIR``:

.. code-block::

  .
  ├── filtered_bam.bam
  ├── filtered_bam_dedup.bam
  ├── filtered_bam_dedup_sorted.bam
  └── filtered_bam_dedup_sorted.bam.bai



Call Sample SNVs
^^^^^^^^^^^^^^^^
Next, you will need to identify SNV genotypes in the pooled bam.

.. code-block:: bash

  singularity exec Demuxafy.sif freebayes -f $FASTA -iXu -C 2 -q 1 $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam > $SCSPLIT_OUTDIR/freebayes_var.vcf
  singularity exec Demuxafy.sif vcftools --gzvcf $SCSPLIT_OUTDIR/freebayes_var.vcf --minQ 30 --recode --recode-INFO-all --out $SCSPLIT_OUTDIR/freebayes_var_qual30

After running these SNV calling steps, you will have the following new files in your ``$SCSPLIT_OUTDIR``:

.. code-block::
  :emphasize-lines: 5,6,7

  .
  ├── filtered_bam.bam
  ├── filtered_bam_dedup.bam
  ├── filtered_bam_dedup_sorted.bam
  ├── filtered_bam_dedup_sorted.bam.bai
  ├── freebayes_var_qual30.log
  ├── freebayes_var_qual30.recode.vcf
  └── freebayes_var.vcf


Demultiplex with ScSplit
^^^^^^^^^^^^^^^^^^^^^^^^
The prepared SNV genotypes and bam file can then be used to demultiplex and call genotypes in each cluster.

.. code-block:: bash

  singularity exec Demuxafy.sif scSplit count -c $VCF -v $SCSPLIT_OUTDIR/freebayes_var_qual30.recode.vcf -i $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam -b $BARCODES -r $SCSPLIT_OUTDIR/ref_filtered.csv -a $SCSPLIT_OUTDIR/alt_filtered.csv -o $SCSPLIT_OUTDIR
  singularity exec Demuxafy.sif scSplit run -r $SCSPLIT_OUTDIR/ref_filtered.csv -a $SCSPLIT_OUTDIR/alt_filtered.csv -n $N -o $SCSPLIT_OUTDIR
  singularity exec Demuxafy.sif scSplit genotype -r $SCSPLIT_OUTDIR/ref_filtered.csv -a $SCSPLIT_OUTDIR/alt_filtered.csv -p $SCSPLIT_OUTDIR/scSplit_P_s_c.csv -o $SCSPLIT_OUTDIR

After running these demultiplexing steps, you will have the following new results:

.. code-block::
  :emphasize-lines: 9,10,11,12,13,14,15,16
  
  .
  ├── alt_filtered.csv
  ├── filtered_bam.bam
  ├── filtered_bam_dedup.bam
  ├── filtered_bam_dedup_sorted.bam
  ├── filtered_bam_dedup_sorted.bam.bai
  ├── freebayes_var_qual30.log
  ├── freebayes_var_qual30.recode.vcf
  ├── freebayes_var.vcf
  ├── ref_filtered.csv
  ├── scSplit_dist_matrix.csv
  ├── scSplit_dist_variants.txt
  ├── scSplit.log
  ├── scSplit_PA_matrix.csv
  ├── scSplit_P_s_c.csv
  ├── scSplit_result.csv
  └── scSplit.vcf

Additional details about outputs are available below in the :ref:`Demuxlet Results and Interpretation <demuxlet-results>`.


ScSplit Summary
^^^^^^^^^^^^^^^
We have provided a script that will provide a summary of the number of droplets classified as doublets, ambiguous and assigned to each cluster by ScSplit_. 
You can run this to get a fast and easy summary of your results.
Just pass the ScSplit_ result file:

.. code-block:: bash

  singularity exec Demuxafy.sif bash scSplit_summary.sh $SCSPLIT_OUTDIR/scSplit_result.csv

which will return the following summary:

  +-----------------+--------------+
  | Classification  | Assignment N |
  +=================+==============+
  | DBL             | 1055         |
  +-----------------+--------------+
  | SNG-0           | 1116         |
  +-----------------+--------------+
  | SNG-10          | 1654         |
  +-----------------+--------------+
  | SNG-11          | 1207         |
  +-----------------+--------------+
  | SNG-12          | 1564         |
  +-----------------+--------------+
  | SNG-13          | 1428         |
  +-----------------+--------------+
  | SNG-14          | 1640         |
  +-----------------+--------------+
  | SNG-2           | 514          |
  +-----------------+--------------+
  | SNG-3           | 1314         |
  +-----------------+--------------+
  | SNG-4           | 1587         |
  +-----------------+--------------+
  | SNG-5           | 1774         |
  +-----------------+--------------+
  | SNG-6           | 1484         |
  +-----------------+--------------+
  | SNG-7           | 1662         |
  +-----------------+--------------+
  | SNG-8           | 1578         |
  +-----------------+--------------+
  | SNG-9           | 1282         |
  +-----------------+--------------+

You can save the summary to file pointing it to the desired output file:

.. code-block:: bash

  singularity exec Demuxafy.sif bash scSplit_summary.sh $SCSPLIT_OUTDIR/scSplit_result.csv > $SCSPLIT_OUTDIR/scSplit_summary.tsv

.. admonition:: Note

  To check if these numbers are consistent with the expected doublet rate in your dataset, you can use our `Doublet Estimation Calculator <test.html>`__.


Correlating Cluster to Donor Reference SNP Genotypes (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have reference SNP genotypes for some or all of the donors in your pool, you can identify which cluster is best correlated with each donor in your reference SNP genotypes. We have provided a script that will do this and provide a heatmap correlation figure and the predicted individual that should be assigned for each cluster. You can either run it with the script by providing the reference SNP genotypes (``$VCF``), the cluster SNP genotypes (``$SCSPLIT_OUTDIR/scSplit.vcf``) and the output directory (``$SCSPLIT_OUTDIR``) You can run this script with:

.. admonition:: Note

  In order to do this, your $VCF must be reference SNP genotypes for the individuals in the pool and cannot be a general vcf with common SNP genotype locations from 1000 Genomes or HRC.

.. tabs::

  .. tab:: With Script

    .. code-block:: bash

      singularity exec Demuxafy.sif Assign_Indiv_by_Geno.R -r $VCF -c $SCSPLIT_OUTDIR/scSplit.vcf -o $SCSPLIT_OUTDIR

    To see the parameter help menu, type:

    .. code-block:: bash

      singularity exec Demuxafy.sif Assign_Indiv_by_Geno.R -h

    Which will print:

    .. code-block:: bash

      usage: Assign_Indiv_by_Geno.R [-h] -r REFERENCE_VCF -c CLUSTER_VCF -o OUTDIR

      optional arguments:
      -h, --help            show this help message and exit
      -r REFERENCE_VCF, --reference_vcf REFERENCE_VCF
                                                      The output directory where results will be saved
      -c CLUSTER_VCF, --cluster_vcf CLUSTER_VCF
                                                      A QC, normalized seurat object with
                                                      classifications/clusters as Idents().
      -o OUTDIR, --outdir OUTDIR
                                                      Number of genes to use in
                                                      'Improved_Seurat_Pre_Process' function.



  .. tab:: Run in R

    You can run the reference vs cluster genotypes manually (possibly because your data doesn't have GT, DS or GP genotype formats) or because you would prefer to alter some of the steps.
    To run the correlations manually, simply start R from the singularity image:

    .. code-block:: R

      singularity exec Demuxafy.sif R

    Once, R has started, you can load the required libraries (included in the singularity image) and run the code.

    .. code-block:: bash

      .libPaths("/usr/local/lib/R/site-library") ### Required so that libraries are loaded from the image instead of locally
      library(tidyr)
      library(tidyverse)
      library(dplyr)
      library(vcfR)
      library(lsa)
      library(ComplexHeatmap)


      ########## Set up paths and variables ##########

      reference_vcf <- "/path/to/reference.vcf"
      cluster_vcf <- "/path/to/scSplit/out/scSplit.vcf"
      outdir <- "/path/to/scSplit/out/"


      ########## Set up functions ##########
      ##### Calculate DS from GP if genotypes in that format #####
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
      ref_geno <- read.vcfR(reference_vcf)
      cluster_geno <- read.vcfR(cluster_vcf)



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
      write_delim(pearson_correlations_out, file = paste0(outdir,"/ref_clust_pearson_correlations.tsv"), delim = "\t" )


      ########## Create correlation figures ##########
      col_fun = colorRampPalette(c("white", "red"))(101)
      pPearsonCorrelations <- Heatmap(as.matrix(pearson_correlations), cluster_rows = T, col = col_fun)

      ########## Save the correlation figures ##########
      png(filename = paste0(outdir,"/ref_clust_pearson_correlation.png"), width = 500)
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

      write_delim(key, file = paste0(outdir,"/Genotype_ID_key.txt"), delim = "\t")



ScSplit Results and Interpretation
----------------------------------
After running the ScSplit_ steps and summarizing the results, you will have a number of files from some of the intermediary steps. Theses are the files that most users will find the most informative:

  - ``scSplit_doublets_singlets.csv``

    - The droplet assignment results. The first column is the droplet barcode and the second column is the droplet type and cluster assignment separated by a dash. For example SNG-9 would indicate that cluster 9 are singlets.

      +--------------------+----------+
      | Barcode            | Cluster  |
      +====================+==========+
      | AAACCTGTCCGAATGT-1 | SNG-0    |
      +--------------------+----------+
      | AAACGGGAGTTGAGAT-1 | SNG-0    |
      +--------------------+----------+
      | AAACGGGCATGTCTCC-1 | SNG-0    |
      +--------------------+----------+
      | AAACGGGTCCACGAAT-1 | SNG-0    |
      +--------------------+----------+
      | AAACGGGTCCAGTAGT-1 | SNG-0    |
      +--------------------+----------+
      | AAACGGGTCGGCTTGG-1 | SNG-0    |
      +--------------------+----------+
      | AAAGATGTCCGAACGC-1 | SNG-0    |
      +--------------------+----------+
      | AAAGATGTCCGTCAAA-1 | SNG-0    |
      +--------------------+----------+
      | AAAGTAGCATCACGTA-1 | SNG-0    |
      +--------------------+----------+
      | ...                | ...      |
      +--------------------+----------+

If you ran the ``Assign_Indiv_by_Geno.R`` script, you will also have the following files:

  - ``Genotype_ID_key.txt``

    - Key of the cluster and assignments for each individual and the Pearson correlation coefficient.

      +-------------+------------+-------------+
      | Genotype_ID | Cluster_ID | Correlation |
      +=============+============+=============+
      | 113_113     | 12         | 0.6448151   |
      +-------------+------------+-------------+
      | 349_350     | 14         | 0.6663323   |
      +-------------+------------+-------------+
      | 352_353     | 7          | 0.6596409   | 
      +-------------+------------+-------------+
      | 39_39       | 6          | 0.6398297   |
      +-------------+------------+-------------+
      | 40_40       | 9          | 0.6191905   |
      +-------------+------------+-------------+
      | 41_41       | 3          | 0.6324396   |
      +-------------+------------+-------------+
      | 42_42       | 4          | 0.6560180   |
      +-------------+------------+-------------+
      | 43_43       | 5          | 0.6672336   |
      +-------------+------------+-------------+
      | 465_466     | 11         | 0.6297396   |
      +-------------+------------+-------------+
      | 596_597     | 13         | 0.6273717   |
      +-------------+------------+-------------+
      | 597_598     | 10         | 0.6627428   |
      +-------------+------------+-------------+
      | 632_633     | 1          | 0.5899685   |
      +-------------+------------+-------------+
      | 633_634     | 0          | 0.6157936   |
      +-------------+------------+-------------+
      | 660_661     | 8          | 0.6423770   |
      +-------------+------------+-------------+

  - ``ref_clust_pearson_correlation.png``

    - Figure of the Pearson correlation coefficients for each cluster-individual pair.

      .. figure:: _figures/OneK1K_scRNA_Sample54_scSplit_pearson_correlation.png

  - ``ref_clust_pearson_correlations.tsv``

    - All of the Pearson correlation coefficients between the clusters and the individuals

      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | Cluster |          113_113    |          349_350    |          352_353    |          39_39      |          40_40      | ... |
      +=========+=====================+=====================+=====================+=====================+=====================+=====+
      | 0       | 0.18419103983986865 | 0.18328230320693129 | 0.19176272973032255 | 0.15376916805897994 | 0.19107524908934623 | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | 1       | 0.19853015287744033 | 0.1981622074955004  | 0.19245840283478327 | 0.17855748333388533 | 0.19455433395443292 | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | 2       | 0.17993959098414505 | 0.15477058833898663 | 0.26412833664924995 | 0.17360648445022142 | 0.16374615160876657 | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | 3       | 0.2128616996153357  | 0.19325148148095284 | 0.21728991668088174 | 0.19346574998787222 | 0.17921651379211084 | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | 4       | 0.17573820413419833 | 0.17629504087312717 | 0.16426156659465307 | 0.17427996983606964 | 0.18322785415879167 | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | ...     | ...                 | ...                 | ...                 | ...                 | ...                 | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+


Merging Results with Other Software Results
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.

Citation
--------
If you used the Demuxafy platform for analysis, please reference our preprint_ as well as `ScSplit <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1852-7>`__.