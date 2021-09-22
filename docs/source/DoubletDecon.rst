.. _DoubletDecon-docs:

DoubletDecon Tutorial
===========================

.. _DoubletDecon: https://github.com/statgen/popscle

DoubletDecon_ is a transcription-based doublet detction software that uses deconvolution to identify doublets using the `R` statistical software.
We have provided a wrapper script that takes common arguments for DoubletDecon_ and also provide example code for you to run manually if you prefer.



Data
----
This is the data that you will need to have preparede to run DoubletDecon_:




Run DoubletDecon
----------------
You can either run DoubletDecon_ with the wrapper script we have provided or you can run it manually if you would prefer to alter more parameters.

.. tabs::

  .. tab:: With Wrapper Script

    .. code-block::

      singularity exec image.sif

  .. tab:: Run in R

    First, you will have to start R.
    We have built R and all the required software to run DoubletDecon_ into the singularity image so you can run it directly from the image.

    .. code-block::

      singularity exec image.sif R

    That will open R in your terminal.
    Next, you can load all the libraries and run DoubletDecon_.

    .. code-block::

      library(DoubletDecon)
      library(tidyverse)
      library(Seurat)
      library(ggplot2)

      ## Set up variables ##
      out <- 
      seurat_object <- "/path/to/preprocessed/seurat_object.rds"




      ## make sure the directory exists ###
      dir.create(out, recursive = TRUE)

      ## Read in Data ##
      seurat <- readRDS(seurat_object)

      ## Preprocess ##
      processed <- Improved_Seurat_Pre_Process(seurat, num_genes=50, write_files=FALSE)

      ## Run Doublet Decon ##
      results <- Main_Doublet_Decon(rawDataFile = processed$newExpressionFile, 
        groupsFile = processed$newGroupsFile, 
        filename = "DoubletDecon_results",
        location = out,
        fullDataFile = NULL, 
        removeCC = FALSE, 
        species = "hsa", 
        rhop = 0.9,
        write = TRUE, 
        PMF = TRUE, 
        useFull = FALSE, 
        heatmap = FALSE, 
        centroids=FALSE, 
        num_doubs=100, 
        only50=FALSE, 
        min_uniq=4, 
        nCores = 1)




      doublets <- fread(paste0(out, "/Final_doublets_groups_DoubletDecon_results.txt"))
      doublets$Barcode <- gsub("\\.", "-",rownames(doublets))
      doublets$DoubletDecon_DropletType <- "doublet"
      doublets$V1 <- NULL
      doublets$V2 <- NULL


      singlets <- fread(paste0(out, "/Final_nondoublets_groups_DoubletDecon_results.txt"))
      singlets$Barcode <- gsub("\\.", "-",rownames(singlets))
      singlets$DoubletDecon_DropletType <- "singlet"
      singlets$V1 <- NULL
      singlets$V2 <- NULL

      doublets_singlets <- rbind(singlets,doublets)

      fwrite(doublets_singlets, paste0(out, "/DoubletDecon_doublets_singlets.tsv", sep = "\t"))


      ### Make a summaruy of the number of singlets and doublets
      summary <- as.data.frame(table(doublets_singlets$DoubletDecon_DropletType))
      colnames(summary) <- c("Classification", "Droplet N")
      write_delim(summary, paste0(out,"/DoubletDecon_doublet_summary.tsv"), "\t")
