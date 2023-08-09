.. _DoubletDecon-docs:

DoubletDecon
===========================

.. _DoubletDecon: https://github.com/EDePasquale/DoubletDecon
.. _preprint: https://www.biorxiv.org/content/10.1101/2022.03.07.483367v1

DoubletDecon_ is a transcription-based doublet detection software that uses deconvolution to identify doublets using the `R` statistical software.
We have provided a wrapper script that takes common arguments for DoubletDecon_ and also provide example code for you to run manually if you prefer.



Data
----
This is the data that you will need to have prepare to run DoubletDecon_:

.. admonition:: Required
  :class: important

  - A QC-filtered and normalized seurat object saved as an ``rds`` object (``$SEURAT_RDS``)

    - For example, using the `Seurat Vignette <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>`__

    - If you run DoubletDecon_ manually, you can use any data format of interest and read in with a method that works for your data.

  - Output directory (``$DOUBLETDECON_OUTDIR``)



Run DoubletDecon
----------------
You can either run DoubletDecon_ with the wrapper script we have provided or you can run it manually if you would prefer to alter more parameters.

.. tabs::

  .. tab:: With Wrapper Script

    .. admonition:: Note

      Since it is hard to predict the correct *rhop* to use for each dataset, we typically run a range.
      For example: 0.6, 0.7, 0.8, 0.9, 1, and 1.1.
      Then we select the results that predict the number of doublets closest to the expected doublet number.
      You can estimate that number with our **doublet calculator**
      The *rhop* parameter can be set with ``-r`` or ``--rhop`` in the command below.

    First, let's assign the variables that will be used to execute each step.

    .. admonition:: Example Variable Settings
      :class: grey

      Below is an example of the variables that we can set up to be used in the command below.
      These are files provided as a :ref:`test dataset <TestData>` available in the :ref:`Data Preparation Documentation <DataPrep-docs>`
      Please replace paths with the full path to data on your system.

      .. code-block:: bash

        DOUBLETDECON_OUTDIR=/path/to/output/DoubletDecon
        SEURAT_RDS=/path/to/TestData_Seurat.rds



    .. code-block:: bash

      singularity exec Demuxafy.sif DoubletDecon.R -o $DOUBLETDECON_OUTDIR -s $SEURAT_RDS

    .. admonition:: HELP! It says my file/directory doesn't exist!
      :class: dropdown

      If you receive an error indicating that a file or directory doesn't exist but you are sure that it does, this is likely an issue arising from Singularity.
      This is easy to fix.
      The issue and solution are explained in detail in the :ref:`Notes About Singularity Images <Singularity-docs>`
      
    You can provide many other parameters as well which can be seen from running a help request:

    .. code-block:: bash

      singularity exec image DoubletDecon.R -h

      usage: DoubletDecon.R [-h] -o OUT -s SEURAT_OBJECT [-g NUM_GENES] [-r RHOP]
                            [-p SPECIES] [-n NCORES] [-c REMOVECC] [-m PMF]
                            [-f HEATMAP] [-t CENTROIDS] [-d NUM_DOUBS] [-5 ONLY50]
                            [-u MIN_UNIQ]

      optional arguments:
        -h, --help            show this help message and exit
        -o OUT, --out OUT     The output directory where results will be saved
        -s SEURAT_OBJECT, --seurat_object SEURAT_OBJECT
                              A QC, normalized seurat object with classifications/clusters as Idents() saved as an rds object.
        -g NUM_GENES, --num_genes NUM_GENES
                              Number of genes to use in 'Improved_Seurat_Pre_Process' function.
        -r RHOP, --rhop RHOP  rhop to use in DoubletDecon - the number of SD from the mean to identify upper limit to blacklist
        -p SPECIES, --species SPECIES
                              The species of your sample. Can be scientific species name, KEGG ID, three letter species abbreviation, or NCBI ID.
        -n NCORES, --nCores NCORES
                              The number of unique cores you would like to use to run DoubletDecon. By default, uses one less than available detected.
        -c REMOVECC, --removeCC REMOVECC
                              Whether to remove clusters enriched in cell cycle genes.
        -m PMF, --pmf PMF     Whether to use unique gene expression in doublet determination.
        -f HEATMAP, --heatmap HEATMAP
                              Whether to generate heatmaps.
        -t CENTROIDS, --centroids CENTROIDS
                              Whether to use centroids instead of medoids for doublet detecting.
        -d NUM_DOUBS, --num_doubs NUM_DOUBS
                              The number of doublets to simulate for each cluster pair.
        -5 ONLY50, --only50 ONLY50
                              Whether to only compute doublets as 50:50 ratio. Default is to use other ratios as well.
        -u MIN_UNIQ, --min_uniq MIN_UNIQ
                              Minimum number of unique genes to rescue a cluster identified as doublets.


  .. tab:: Run in R

    First, you will have to start R.
    We have built R and all the required software to run DoubletDecon_ into the singularity image so you can run it directly from the image.

    .. code-block:: bash

      singularity exec Demuxafy.sif R

    That will open R in your terminal.
    Next, you can load all the libraries and run DoubletDecon_.

    .. code-block:: R

      .libPaths("/usr/local/lib/R/site-library") ### This is required so that R uses the libraries loaded in the image and not any local libraries
      library(DoubletDecon)
      library(tidyverse)
      library(Seurat)
      library(ggplot2)
      library(data.table)

      ## Set up variables ##
      out <- "/path/to/doubletdecon/outdir"
      SEURAT_RDSect <- "/path/to/preprocessed/SEURAT_RDSect.rds"




      ## make sure the directory exists ###
      dir.create(out, recursive = TRUE)

      ## Read in Data ##
      seurat <- readRDS(SEURAT_RDSect)

      ## Preprocess ##
      processed <- Improved_Seurat_Pre_Process(seurat, num_genes=50, write_files=FALSE)

      ## Run Doublet Decon ##
      results <- Main_Doublet_Decon(rawDataFile = processed$newExpressionFile, 
        groupsFile = processed$newGroupsFile, 
        filename = "DoubletDecon_results",
        location = paste0(out, "/"),
        fullDataFile = NULL, 
        removeCC = FALSE, 
        species = "hsa", 
        rhop = 0.9,                         ## We recommend testing multiple rhop parameters to find which fits your data the best
        write = TRUE, 
        PMF = TRUE, 
        useFull = FALSE,
        heatmap = FALSE, 
        centroids=FALSE, 
        num_doubs=100, 
        only50=FALSE, 
        min_uniq=4, 
        nCores = 1)




      doublets <- read.table(paste0(out, "/Final_doublets_groups_DoubletDecon_results.txt"))
      doublets$Barcode <- gsub("\\.", "-",rownames(doublets))
      doublets$DoubletDecon_DropletType <- "doublet"
      doublets$V1 <- NULL
      doublets$V2 <- NULL


      singlets <- read.table(paste0(out, "/Final_nondoublets_groups_DoubletDecon_results.txt"))
      singlets$Barcode <- gsub("\\.", "-",rownames(singlets))
      singlets$DoubletDecon_DropletType <- "singlet"
      singlets$V1 <- NULL
      singlets$V2 <- NULL

      doublets_singlets <- rbind(singlets,doublets)

      fwrite(doublets_singlets, paste0(out, "/DoubletDecon_doublets_singlets.tsv"), sep = "\t", append = FALSE)


      ### Make a summary of the number of singlets and doublets
      summary <- as.data.frame(table(doublets_singlets$DoubletDecon_DropletType))
      colnames(summary) <- c("Classification", "Droplet N")
      fwrite(summary, paste0(out,"/DoubletDecon_doublet_summary.tsv"), sep = "\t", append = FALSE)



DoubletDecon Results and Interpretation
----------------------------------------
After running the DoubletDecon_, you will have multiple files in the ``$DOUBLETDECON_OUTDIR``:

.. code-block:: bash

  /path/to/output/DoubletDecon
  ├── data_processed_DoubletDecon_results.txt
  ├── data_processed_reclust_DoubletDecon_results.txt
  ├── DoubletDecon_doublets_singlets.tsv
  ├── DoubletDecon_doublet_summary.tsv
  ├── DoubletDecon_results.log
  ├── DRS_doublet_table_DoubletDecon_results.txt
  ├── DRS_results_DoubletDecon_results.txt
  ├── Final_doublets_exp_DoubletDecon_results.txt
  ├── Final_doublets_groups_DoubletDecon_results.txt
  ├── Final_nondoublets_exp_DoubletDecon_results.txt
  ├── Final_nondoublets_groups_DoubletDecon_results.txt
  ├── groups_processed_DoubletDecon_results.txt
  ├── groups_processed_reclust_DoubletDecon_results.txt
  ├── new_PMF_results_DoubletDecon_results.txt
  ├── resultsreadable_synths.txt
  └── Synth_doublet_info_DoubletDecon_results.txt


DoubletDecon_ puts most of the results in multiple separate files. 
However, the wrapper script and the example code has some steps to combine these results together into a single file, which will likely be the most informative output.
These are the files that we think will be the most helpful for users:

- ``DoubletDecon_doublet_summary.tsv``
  
  - A summary of the number of singlets and doublets predicted by DoubletDecon_.

    +----------------+-----------+
    |Classification  | Droplet N |
    +================+===========+
    |doublet         | 1510      |
    +----------------+-----------+
    |singlet         | 19470     |
    +----------------+-----------+

    - To check whether the numbe of doublets identified by DoubletDecon_ is consistent with the expected doublet rate expected based on the number of droplets that you captured, you can use our `Expected Doublet Estimation Calculator <test.html>`__.

- ``DoubletDecon_doublets_singlets.tsv``

  - The per-barcode singlet and doublet classification from DoubletDecon_.

    +-------------------------+--------------------------+
    | Barcode                 | DoubletDecon_DropletType |
    +=========================+==========================+
    | AAACCTGAGCAGCGTA-1      | singlet                  |
    +-------------------------+--------------------------+
    | AAACCTGAGCGATGAC-1      | singlet                  |
    +-------------------------+--------------------------+
    | AAACCTGAGCGTAGTG-1      | singlet                  |
    +-------------------------+--------------------------+
    | AAACCTGAGGCTCATT-1      | singlet                  |
    +-------------------------+--------------------------+
    | AAACCTGAGTAGCCGA-1      | singlet                  |
    +-------------------------+--------------------------+
    | ...                     | ...                      |
    +-------------------------+--------------------------+


Merging Results with Other Software Results
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.


Citation
--------
If you used the Demuxafy platform for analysis, please reference our preprint_ as well as `DoubletDecon <https://www.sciencedirect.com/science/article/pii/S2211124719312860>`__.