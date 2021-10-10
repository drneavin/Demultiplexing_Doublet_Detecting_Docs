.. _scDblFinder-docs:

scDblFinder Tutorial
===========================

.. _scDblFinder: https://github.com/plger/scDblFinder

scDblFinder is a transcriptome-based doublet detecting method that uses doublet simulation from droplets in the dataset to identify doublets.
We have provided a wrapper script that takes common arguments for scDblFinder_ and also provide example code for you to run manually if you prefer.



Data
----
This is the data that you will need to have preparede to run scDblFinder_:

.. admonition:: Required
  :class: important

  - A counts matrix (``$MATRIX``)
  
    - DoubletDetection expects counts to be in the cellranger output format (directory containint ``barcodes.tsv``, ``genes.tsv`` and ``matrix.mtx`` **or** ``barcodes.tsv.gz``, ``features.tsv.gz`` and ``matrix.mtx.gz``)

	  - If you don't have your data in this format, you can run scDblFinder_ manually in python and load the data in using a method of your choosing.

  - Output directory (``$SCDBLFINDER_OUTDIR``)

    - If you don't provide an ``$SCDBLFINDER_OUTDIR``, the results will be written to the present working directory.




Run scDblFinder
----------------
You can either run scDblFinder_ with the wrapper script we have provided or you can run it manually if you would prefer to alter more parameters.

.. tabs::

  .. tab:: With Wrapper Script

    .. code-block:: bash

		  singularity exec Demuxafy.sif Rscript scDblFinder.R -o $SCDBLFINDER_OUTDIR -t $MATRIX


  .. tab:: Run in R

    First, you will have to start R.
    We have built R and all the required software to run scDblFinder_ into the singularity image so you can run it directly from the image.

    .. code-block:: bash

      singularity exec Demuxafy.sif R


    That will open R in your terminal.
    Next, you can load all the libraries and run scDblFinder_.

    .. code-block:: R

      .libPaths("/usr/local/lib/R/site-library") ### This is required so that R uses the libraries loaded in the image and not any local libraries
      library(scDblFinder)
      library(Seurat)
      library(SingleCellExperiment)
      library(tidyverse)

      ## Set up variables and parameters ##
      out <- "/path/to/scds/outdir/"
      tenX_matrix <- "/path/to/counts/matrix/dir/"

      dir.create(out, recursive = TRUE)
      print(paste0("Using the following counts directory: ", tenX_matrix))



      ### Read in data as an sce object ###
      counts <- Read10X(tenX_matrix, gene.column = 1)
      sce <- SingleCellExperiment(list(counts=counts))


      ## Calculate doublet ratio ###
      doublet_ratio <- ncol(sce)/1000*0.008


      ### Calculate Singlets and Doublets ###
      sce <- scDblFinder(sce, dbr=doublet_ratio)


      
      ### Make a dataframe of the results ###
      results <- data.frame("Barcode" = rownames(colData(sce)), "scDblFinder_DropletType" = sce$scDblFinder.class, "scDblFinder_Score" = sce$scDblFinder.score)


      write_delim(results, path = paste0(out,"/scDblFinder_doublets_singlets.tsv"), delim = "\t")

      ### Calculate number of doublets and singlets ###
      summary <- as.data.frame(table(results$scDblFinder_DropletType))
      colnames(summary) <- c("Classification", "Droplet N")
      write_delim(summary, paste0(out,"/scDblFinder_doublet_summary.tsv"), "\t")



scDblFinder Results and Interpretation
----------------------------------------
After running the scDblFinder_ with the wrapper script or manually you should have two files in the ``$SCDBLFINDER_OUTDIR``:

.. code-block:: bash

	.
	├── scDblFinder_doublets_singlets.tsv
	└── scDblFinder_doublet_summary.tsv

Here's a more detaild description of each of those files:

- ``scDblFinder_doublet_summary.tsv``

  - A sumamry of the number of singlets and doublets predicted by scDblFinder_.

    +----------------+-----------+
    |Classification  | Droplet N |
    +================+===========+
    |doublet         | 3323      |
    +----------------+-----------+
    |singlet         | 17659     |
    +----------------+-----------+

    - To check whether the numbe of doublets identified by scDblFinder_ is consistent with the expected doublet rate expected based on the number of droplets that you captured, you can use our `Expected Doublet Estimation Calculator <test.html>`__.

- ``scDblFinder_doublets_singlets.tsv``

  - The per-barcode singlet and doublet classification from scDblFinder_.

    +-------------------------+-------------------------+--------------------------+
    | Barcode                 | scDblFinder_DropletType | scDblFinder_Score        |
    +=========================+=========================+==========================+
    | AAACCTGAGATAGCAT-1      | singlet                 | 0.0033526041079312563    |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGCAGCGTA-1      | doublet                 | 0.9937564134597778       |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGCGATGAC-1      | singlet                 | 5.045032594352961e-      |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGCGTAGTG-1      | singlet                 | 0.007504515815526247     |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGGAGTTTA-1      | singlet                 | 0.00835108570754528      |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGGCTCATT-1      | singlet                 | 0.028838597238063812     |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGGGCACTA-1      | doublet                 | 0.9985504746437073       |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGTAATCCC-1      | singlet                 | 0.005869860760867596     |
    +-------------------------+-------------------------+--------------------------+
    | ...                     | ...                     | ...                      |
    +-------------------------+-------------------------+--------------------------+


Merging Results with Other Software Restults
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.

Citation
--------
If you used this workflow for analysis, please reference our paper (REFERENCE) as well as `scDblFinder <https://github.com/plger/scDblFinder>`__.

