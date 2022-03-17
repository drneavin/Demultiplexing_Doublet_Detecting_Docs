.. _Souporcell-docs:

Souporcell
===========================

.. _Souporcell: https://github.com/wheaton5/souporcell
.. _preprint: https://www.biorxiv.org/content/10.1101/2022.03.07.483367v1

Souporcell_ is a genotype-free demultiplexing software that does not require you to have SNP genotypes the donors in your multiplexed capture.
However, it can natively integrate SNP genotypes into the demultiplexing if you have them available for **all** the donors in your pool.
If you don't have the reference SNP genotypes for all the donors in your multiplexed pool, we have provided some scripts that will help identify clusters from given donors after running Souporcell_ without the SNP genotypes.
Depending on your downstream analyses, if you have reference SNP genotypes for donors in your pool, you could also use :ref:`Demuxlet <Demuxlet-docs>`, or :ref:`Vireo<Vireo-docs>`.

One advantage that we have found immensely helpful about Souporcell_ is that it provides an ambient RNA estimate for the pool.
This can be helpful to identify samples that may have high ambient RNA estimates early in the analysis pipeline so that it can be accounted for throughout downstream analyses.






Data
----
This is the data that you will need to have prepare to run Souporcell_:

.. admonition:: Required
  :class: important

  - Common SNP genotypes vcf (``$VCF``)

    - While not exactly required, using common SNP genotype locations enhances accuracy

      - If you have reference SNP genotypes for individuals in your pool, you can use those

      - If you do not have reference SNP genotypes, they can be from any large population resource (i.e. 1000 Genomes or HRC)

    - Filter for common SNPs (> 5% minor allele frequency) and SNPs overlapping genes

  - Barcode file (``$BARCODES``)

  - Number of samples in pool (``$N``)
  
  - Bam file (``$BAM``)

    - Aligned single cell reads

  - Reference fasta (``$FASTA``)
  
    - that your reads were aligned to (or at least the same genome)

  - Output directory (``$SOUPORCELL_OUTDIR``)



Run Souporcell
--------------
You can run Souporcell_ with or without reference SNP genotypes - follow the instructions for each bellow:


.. tabs::

  .. tab:: Without Reference SNP Genotypes

    If you don't have reference SNP genotypes for all of your donors, you can run souporcell with the following command, providing an appropriate thread number (``$THREADS``) for your system .
    Don't worry if you only have reference SNP genotypes for a subset of your donors, we have a script that will correlate the cluster and reference SNP genotypes.

    .. code-block:: bash

      singularity exec Demuxafy.sif souporcell_pipeline.py -i $BAM -b $BARCODES -f $FASTA -t $THREADS -o $SOUPORCELL_OUTDIR -k $N --common_variants $VCF

  .. tab:: With Reference SNP Genotypes

    If you have reference SNP genotypes for **all** of your donors, you can run souporcell with the following command, providing an appropriate thread number (``$THREADS``) for your system and listing the donor ids that correspond in the ``$VCF`` file

    .. code-block:: bash

      singularity exec Demuxafy.sif souporcell_pipeline.py -i $BAM -b $BARCODES -f $FASTA -t $THREADS -o $SOUPORCELL_OUTDIR -k $N --known_genotypes $VCF --known_genotypes_sample_names donor1 donor donor3 donor4

    .. admonition:: Note

      Souporcell can currently only be executed when either **all** or **none** of the individuals that have been pooled have SNP genotypes.
      Further, the output still has cluster numbers but they should correspond to the order that you listed your individuals.
      For example, if you have two individuals in your pool (donorA and donorB) and input them as ``--known_genotypes_sample_names donorA donorB``, then the output will have two clusters: 0 and 1.
      donorA will correspond to 0 and donorB will correspond to 1.

      Even when we have reference SNP genotypes, we typically runn Souporcell without reference SNP genotypes and then use the cluster vs individual correlations (below) to assign clusters to individuals.

If Souporcell_ is successful, you will have these files in your ``$SOUPORCELL_OUTDIR``:

.. code-block:: bash

  .
  ├── alt.mtx
  ├── ambient_rna.txt
  ├── cluster_genotypes.vcf
  ├── clustering.done
  ├── clusters.err
  ├── clusters_tmp.tsv
  ├── clusters.tsv
  ├── common_variants_covered_tmp.vcf
  ├── common_variants_covered.vcf
  ├── consensus.done
  ├── depth_merged.bed
  ├── doublets.err
  ├── fastqs.done
  ├── minimap.err
  ├── ref.mtx
  ├── remapping.done
  ├── retag.err
  ├── retagging.done
  ├── souporcell_minimap_tagged_sorted.bam
  ├── souporcell_minimap_tagged_sorted.bam.bai
  ├── troublet.done
  ├── variants.done
  └── vartrix.done

Additional details about outputs are available below in the :ref:`Souporcell Results and Interpretation <souporcell-results>`.



Souporcell Summary
^^^^^^^^^^^^^^^^^^
We have provided a script that will provide a summary of the number of droplets classified as doublets, ambiguous and assigned to each cluster by Souporcell_. 
You can run this to get a fast and easy summary of your results by providing the souporcell result file:

.. code-block:: bash

  singularity exec Demuxafy.sif bash souporcell_summary.sh $SOUPORCELL_OUTDIR/clusters.tsv

which should print:

  +-----------------+--------------+
  | Classification  | Assignment N |
  +=================+==============+
  | 0               | 1441         |
  +-----------------+--------------+
  | 1               | 980          |
  +-----------------+--------------+
  | 10              | 1285         |
  +-----------------+--------------+
  | 11              | 1107         |
  +-----------------+--------------+
  | 12              | 1315         |
  +-----------------+--------------+
  | 13              | 1529         |
  +-----------------+--------------+
  | 2               | 1629         |
  +-----------------+--------------+
  | 3               | 1473         |
  +-----------------+--------------+
  | 4               | 1381         |
  +-----------------+--------------+
  | 5               | 1360         |
  +-----------------+--------------+
  | 6               | 1157         |
  +-----------------+--------------+
  | 7               | 892          |
  +-----------------+--------------+
  | 8               | 1111         |
  +-----------------+--------------+
  | 9               | 1565         |
  +-----------------+--------------+
  | doublet         | 2757         |
  +-----------------+--------------+

or you can write the results to file:

.. code-block::

  singularity exec Demuxafy.sif bash souporcell_summary.sh $SOUPORCELL_OUTDIR/clusters.tsv > $SOUPORCELL_OUTDIR/souporcell_summary.tsv


.. admonition:: Note

  To check if these numbers are consistent with the expected doublet rate in your dataset, you can use our `Doublet Estimation Calculator <test.html>`__.




If the souporcell summary is successful, you will have this new file in your ``$SOUPORCELL_OUTDIR``:

.. code-block:: bash
  :emphasize-lines: 21

  .
  ├── alt.mtx
  ├── ambient_rna.txt
  ├── cluster_genotypes.vcf
  ├── clustering.done
  ├── clusters.err
  ├── clusters_tmp.tsv
  ├── clusters.tsv
  ├── common_variants_covered_tmp.vcf
  ├── common_variants_covered.vcf
  ├── consensus.done
  ├── depth_merged.bed
  ├── doublets.err
  ├── fastqs.done
  ├── minimap.err
  ├── ref.mtx
  ├── remapping.done
  ├── retag.err
  ├── retagging.done
  ├── souporcell_minimap_tagged_sorted.bam
  ├── souporcell_summary.tsv
  ├── troublet.done
  ├── variants.done
  └── vartrix.done

Additional details about outputs are available below in the :ref:`Souporcell Results and Interpretation <souporcell-results>`.


.. _souporcell-results:

Correlating Cluster to Donor Reference SNP Genotypes (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have reference SNP genotypes for some or all of the donors in your pool, you can identify which cluster is best correlated with each donor in your reference SNP genotypes. We have provided a script that will do this and provide a heatmap correlation figure and the predicted individual that should be assigned for each cluster. You can either run it with the script by providing the reference SNP genotypes (``$VCF``), the cluster SNP genotypes (``$SOUPORCELL_OUTDIR/cluster_genotypes.vcf``) and the output directory (``$SOUPORCELL_OUTDIR``) You can run this script with:

.. admonition:: Note

  In order to do this, your $VCF must be reference SNP genotypes for the individuals in the pool and cannot be a general vcf with common SNP genotype locations from 1000 Genomes or HRC.

.. tabs::

  .. tab:: With Script

    .. code-block:: bash

      singularity exec Demuxafy.sif Assign_Indiv_by_Geno.R -r $VCF -c $SOUPORCELL_OUTDIR/cluster_genotypes.vcf -o $SOUPORCELL_OUTDIR

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
      cluster_vcf <- "/path/to/souporcell/out/cluster_genotypes.vcf"
      outdir <- "/path/to/souporcell/out/"


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


After correlating the cluster and reference donor SNP genotypes, you should have the new results in your directory:


If the souporcell summary is successful, you will have this new file in your ``$SOUPORCELL_OUTDIR``:

.. code-block:: bash
  :emphasize-lines: 15,16,19,20

  .
  ├── alt.mtx
  ├── ambient_rna.txt
  ├── cluster_genotypes.vcf
  ├── clustering.done
  ├── clusters.err
  ├── clusters_tmp.tsv
  ├── clusters.tsv
  ├── common_variants_covered_tmp.vcf
  ├── common_variants_covered.vcf
  ├── consensus.done
  ├── depth_merged.bed
  ├── doublets.err
  ├── fastqs.done
  ├── Genotype_ID_key.txt
  ├── Individual_genotypes_subset.vcf.gz
  ├── minimap.err
  ├── ref.mtx
  ├── ref_clust_pearson_correlation.png
  ├── ref_clust_pearson_correlations.tsv
  ├── remapping.done
  ├── retag.err
  ├── retagging.done
  ├── souporcell_minimap_tagged_sorted.bam
  ├── souporcell_summary.tsv
  ├── troublet.done
  ├── variants.done
  └── vartrix.done

Additional details about outputs are available below in the :ref:`Souporcell Results and Interpretation <souporcell-results>`.



Souporcell Results and Interpretation
-------------------------------------
After running the Souporcell_ steps and summarizing the results, you will have a number of files from some of the intermediary steps. 
These are the files that most users will find the most informative:


    - To check if these numbers are consistent with the expected doublet rate in your dataset, you can use our `Expected Doublet Estimation Calculator <test.html>`__.

  - ``clusters.tsv``

    - The Souporcell_ droplet classifications with the log probabilities of each donor and doublet vs singlet.

      +-------------------------+----------+-----------------+-------------------------+-------------------------+-------------------------+--------------------------+-------------------------+-------------------------+----------------------+-------------------------+---------------------------+-------------------------+------------------------+-------------------------+-------------------------+--------------------------+------------------------+----------------------+
      | barcode                 | status   | assignment      | log_prob_singleton      | log_prob_doublet        | cluster0                |cluster1                  | cluster2                | cluster3                | cluster4             | cluster5                | cluster6                  | cluster7                | cluster8               |cluster9                 | cluster10               | cluster11                | cluster12              | cluster13            |
      +-------------------------+----------+-----------------+-------------------------+-------------------------+-------------------------+--------------------------+-------------------------+-------------------------+----------------------+-------------------------+---------------------------+-------------------------+------------------------+-------------------------+-------------------------+--------------------------+------------------------+----------------------+
      | AAACCTGAGATAGCAT-1      | singlet  | 6               | -47.4906809612613       | -67.16353115825044      | -189.38489711217204     |-167.22863078578243       | -175.6243866125455      | -195.88836978493757     | -147.1278571646738   | -162.71464140958287     | -47.4906809612613         | -147.57558470556503     |-142.24543450475267     | -137.94217556189426     | -171.6924681433834      | -192.9070590872178       |-162.2042834302814      | -141.9657291979218   |
      +-------------------------+----------+-----------------+-------------------------+-------------------------+-------------------------+--------------------------+-------------------------+-------------------------+----------------------+-------------------------+---------------------------+-------------------------+------------------------+-------------------------+-------------------------+--------------------------+------------------------+----------------------+
      | AAACCTGAGCAGCGTA-1      | singlet  | 11              | -102.80051804401324     | -158.38006105671326     | -357.5113573904763      |-403.04676141772245       | -465.3312627534814      | -368.72445203224066     | -362.5022337777086   | -377.5322002577741      | -400.12257643517944       | -436.7935123280712      |-364.36305907429954     | -434.8878131790703      | -393.42953156344277     | -102.80051804401324      |-369.5775718688619      | -403.83637627549155  |
      +-------------------------+----------+-----------------+-------------------------+-------------------------+-------------------------+--------------------------+-------------------------+-------------------------+----------------------+-------------------------+---------------------------+-------------------------+------------------------+-------------------------+-------------------------+--------------------------+------------------------+----------------------+
      | AAACCTGAGCGATGAC-1      | singlet  | 5               | -39.97694257579923      | -53.76617956926222      | -135.58935896223636     |-129.29863536547518       | -122.20920829636167     | -99.54420652897485      | -139.8403265674046   | -39.97694257579923      | -136.5313839118704        | -139.57805752070823     |-113.63185227373309     | -117.89083888468238     | -126.95555633151154     | -167.2476854256994       |-127.05455963457722     | -123.63808626520557  |
      +-------------------------+----------+-----------------+-------------------------+-------------------------+-------------------------+--------------------------+-------------------------+-------------------------+----------------------+-------------------------+---------------------------+-------------------------+------------------------+-------------------------+-------------------------+--------------------------+------------------------+----------------------+
      | AAACCTGAGCGTAGTG-1      | singlet  | 3               | -66.73447359908208      | -79.59130566934348      | -146.47954690347862     |-197.54291944344263       | -211.47148694945332     | -66.73447359908208      | -163.94180016636983  | -173.4754549428176      | -183.73592914945144       | -163.7126225130574      |-172.5171380662907      | -231.65011940831332     | -197.42816500995383     | -167.68988627905136      |-165.7006532267023      | -174.74052654720117  |
      +-------------------------+----------+-----------------+-------------------------+-------------------------+-------------------------+--------------------------+-------------------------+-------------------------+----------------------+-------------------------+---------------------------+-------------------------+------------------------+-------------------------+-------------------------+--------------------------+------------------------+----------------------+
      | ...                     | ...      |                 | ...                     | ...                     | ...                     |...                       | ...                     | ...                     | ...                  | ...                     | ...                       | ...                     |...                     | ...                     | ...                     | ...                      |...                     | ...                  |
      +-------------------------+----------+-----------------+-------------------------+-------------------------+-------------------------+--------------------------+-------------------------+-------------------------+----------------------+-------------------------+---------------------------+-------------------------+------------------------+-------------------------+-------------------------+--------------------------+------------------------+----------------------+
      

  - ``ambient_rna.txt``

    - The estimated ambient RNA percent in the pool. We typically see < 5% for scRNA-seq PBMCs and < 10% for other scRNA-seq cell types.

      .. code-block:: bash

        ambient RNA estimated as 4.071468697320357%

  
If you ran the ``Assign_Indiv_by_Geno.R`` script, you will also have the following files:

  - ``Genotype_ID_key.txt``

    - Key of the cluster and assignments for each individual and the Pearson correlation coefficient.

      +-------------+------------+-------------+
      | Genotype_ID | Cluster_ID | Correlation |
      +=============+============+=============+
      | 113_113     |  5         | 0.9365902   |
      +-------------+------------+-------------+
      | 349_350     |  3         | 0.9484794   |
      +-------------+------------+-------------+
      | 352_353     |  2         | 0.9385500   | 
      +-------------+------------+-------------+
      | 39_39       |  12        | 0.9325007   |
      +-------------+------------+-------------+
      | 40_40       |  8         | 0.9252865   |
      +-------------+------------+-------------+
      | 41_41       |  6         | 0.9282633   |
      +-------------+------------+-------------+
      | 42_42       |  0         | 0.9387788   |
      +-------------+------------+-------------+
      | 43_43       |  9         | 0.9497327   |
      +-------------+------------+-------------+
      | 465_466     |  11        | 0.9234109   |
      +-------------+------------+-------------+
      | 596_597     |  10        | 0.9277824   |
      +-------------+------------+-------------+
      | 597_598     |  13        | 0.9435752   |
      +-------------+------------+-------------+
      | 632_633     |  7         | 0.9179054   |
      +-------------+------------+-------------+
      | 633_634     |  1         | 0.9222734   |
      +-------------+------------+-------------+
      | 660_661     |  4         | 0.9368751   |
      +-------------+------------+-------------+


  - ``ref_clust_pearson_correlation.png``

    - Figure of the Pearson correlation coefficients for each cluster-individual pair.

      .. figure:: _figures/OneK1K_scRNA_Sample54_souporcell_pearson_correlation.png

  - ``ref_clust_pearson_correlations.tsv``

    - All of the Pearson correlation coefficients between the clusters and the individuals

      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | Cluster |          113_113    |          349_350    |          352_353    |          39_39      |          40_40      | ... |
      +=========+=====================+=====================+=====================+=====================+=====================+=====+
      | 0       | 0.4578087241392215  | 0.4589573335017816  | 0.46351292453350446 | 0.48926720614880104 | 0.4841871441103791  | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | 1       | 0.45706434043842825 | 0.48280445273461425 | 0.4702618797322548  | 0.4678187806965093  | 0.4801164797099736  | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | 2       | 0.4760176832308062  | 0.45281488606508186 | 0.9385500036660724  | 0.47703829279476667 | 0.47639771569917855 | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | 3       | 0.4771709808299328  | 0.9484794352067363  | 0.4598361363766827  | 0.4698832593827229  | 0.4822779587579728  | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | 4       | 0.4851872933346752  | 0.48480637867431775 | 0.4908275654324142  | 0.48900594491809124 | 0.4647100675599844  | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+
      | ...     | ...                 | ...                 | ...                 | ...                 | ...                 | ... |
      +---------+---------------------+---------------------+---------------------+---------------------+---------------------+-----+


Merging Results with Other Software Results
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.


Citation
--------
If you used the Demuxafy platform for analysis, please reference our preprint_ as well as `Souporcell <https://www.nature.com/articles/s41592-020-0820-1>`__.