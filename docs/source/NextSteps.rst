.. _NextSteps-docs:

Next Steps
==============

.. _publication: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03224-8

After you've run the Demultiplexing and/or Doublet Detecting softwares you would like, you can easily add the results to your single cell data with your analysis program of choice.
This works best using the `combined results <Combine-docs>`__ files since they are well formatted for this function.
Some common softwares used to analyze single cell data include `Seurat <https://satijalab.org/seurat/>`__ (in R), `Scanpy <https://scanpy.readthedocs.io/en/stable/>`__ (in python) and the `Loupe Browser <https://www.10xgenomics.com/products/loupe-browser>`__ (`10x Genomics <https://www.10xgenomics.com/>`__).

We've provided some example methods to add these results to single cell data structures in with each of these packages below.


Seurat
---------

The results can be added to a `Seurat <https://satijalab.org/seurat/>`__ object in R.

First, open R from the Singularity image:


.. code-block:: bash

  singularity exec Demuxafy.sif R

This is some basic code in R that will add the combined results to a `Seurat <https://satijalab.org/seurat/>`__ object.

.. code-block:: R

  .libPaths("/usr/local/lib/R/site-library") ### This is required so that R uses the libraries loaded in the image and not any local libraries
  library(Seurat)
  library(tidyverse)

  ## Read in the data
  counts_matrix <- Read10X(data.dir = "/path/to/10x/matrix/directory/")

  ## Create a seurat object that contains the counts
  seurat <- CreateSeuratObject(counts = counts_matrix, min.cells = 3, min.features = 200)

  ## Read in the demultiplexing and doublet detecting results
  demuxafy <- read.table("/path/to/combined/results/combined_results_w_combined_assignments.tsv", sep = "\t", header=TRUE)
  rownames(demuxafy) <- demuxafy$Barcode

  ## Add the demuxafy data to the Seurat object
  seurat <- AddMetaData(seurat, demuxafy)

  ## Check that the data was correctly added
  head(seurat@meta.data)


If the data was correctly added to the Seurat object, you should be able to see the data in the meta.data slot:

.. code-block:: R

                        orig.ident nCount_RNA nFeature_RNA            Barcode
  CGTTAGATCTAGAGTC-1 SeuratProject       3001          934 CGTTAGATCTAGAGTC-1
  CAGGTGCAGGTCATCT-1 SeuratProject       1778          658 CAGGTGCAGGTCATCT-1
  GTGCATAGTATAGGGC-1 SeuratProject       2224          925 GTGCATAGTATAGGGC-1
  TTCGAAGTCCAGTAGT-1 SeuratProject       1268          428 TTCGAAGTCCAGTAGT-1
  TAGTTGGGTCTCATCC-1 SeuratProject        981          476 TAGTTGGGTCTCATCC-1
  CGCGGTAAGGATTCGG-1 SeuratProject       2982          710 CGCGGTAAGGATTCGG-1
                    Souporcell_DropletType Souporcell_Cluster
  CGTTAGATCTAGAGTC-1                singlet                  7
  CAGGTGCAGGTCATCT-1                singlet                  7
  GTGCATAGTATAGGGC-1                singlet                  7
  TTCGAAGTCCAGTAGT-1                singlet                  7
  TAGTTGGGTCTCATCC-1                singlet                  7
  CGCGGTAAGGATTCGG-1                singlet                  7
                    Souporcell_Individual_Assignment scds_score scds_DropletType
  CGTTAGATCTAGAGTC-1                                7 0.70178412          singlet
  CAGGTGCAGGTCATCT-1                                7 0.04977119          singlet
  GTGCATAGTATAGGGC-1                                7 0.26765665          singlet
  TTCGAAGTCCAGTAGT-1                                7 0.06354318          singlet
  TAGTTGGGTCTCATCC-1                                7 0.12028167          singlet
  CGCGGTAAGGATTCGG-1                                7 0.11765345          singlet
                    solo_DropletType solo_DropletScore
  CGTTAGATCTAGAGTC-1          singlet           doublet
  CAGGTGCAGGTCATCT-1          singlet           doublet
  GTGCATAGTATAGGGC-1          singlet           doublet
  TTCGAAGTCCAGTAGT-1          singlet           doublet
  TAGTTGGGTCTCATCC-1          singlet           doublet
  CGCGGTAAGGATTCGG-1          singlet           doublet
                    MajoritySinglet_DropletType
  CGTTAGATCTAGAGTC-1                     singlet
  CAGGTGCAGGTCATCT-1                     singlet
  GTGCATAGTATAGGGC-1                     singlet
  TTCGAAGTCCAGTAGT-1                     singlet
  TAGTTGGGTCTCATCC-1                     singlet
  CGCGGTAAGGATTCGG-1                     singlet
                    MajoritySinglet_Individual_Assignment
  CGTTAGATCTAGAGTC-1                                     7
  CAGGTGCAGGTCATCT-1                                     7
  GTGCATAGTATAGGGC-1                                     7
  TTCGAAGTCCAGTAGT-1                                     7
  TAGTTGGGTCTCATCC-1                                     7
  CGCGGTAAGGATTCGG-1                                     7




Scanpy
-----------

The results can be added to a AnnData object for analysis with `Scanpy <https://scanpy.readthedocs.io/en/stable/>`__.

First, open python from the Singularity image:


.. code-block:: bash

  singularity exec Demuxafy.sif python

This is some basic code in python that will add the combined results to a `Scanpy <https://scanpy.readthedocs.io/en/stable/>`__.

.. code-block:: python

  import pandas as pd
  import scanpy as sc
  import numpy as np

  ### Read in the data to an AnnData object
  adata = sc.read_10x_mtx("/path/to/10x/matrix/directory/")

  ### Read in the demultiplexing and doublet detecting results
  demuxafy = pd.read_table("/path/to/combined/results/combined_results.tsv", sep="\t")

  ### Filter the AnnData object for droplet barcodes
  adata = adata[np.isin(adata.obs.index,demuxafy["Barcode"])]

  ### Order the demuxafy droplets in same orderas the AnnDataa
  adata_obs = pd.DataFrame(adata.obs)
  adata_obs['Barcode'] = adata.obs.index

  demuxafy_ordered = adata_obs.merge(demuxafy, on = "Barcode")
  demuxafy_ordered.index = demuxafy_ordered["Barcode"]

  ### Add demuxafy data to the AnnData
  adata.obs = demuxafy_ordered


Loupe
--------------

The Demuxafy results from Combine_Results.R can be directly uploaded to the Loupe browser in the 'Categories mode'.
Simpley select 'Import Categories' and select the ``combined_results.tsv`` file to upload it and explore the annotation on your data.
More detailed instructions are provided by `10x Genomics <https://www.10xgenomics.com/>`__ in the 'Categories mode' section of their `Software Support <https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/tutorial-navigation#view-selector>`__




Citation
--------
If you used the Demuxafy platform for analysis, please reference our publication_.