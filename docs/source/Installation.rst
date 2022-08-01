Installation
==========================
Installation should be pretty painless (we hope).
We have  provided all the softwares in a singularity image which provides continuity across different computing platforms (see `HPCNG Singluarity <https://singularity.hpcng.org/>`__ and `Sylabs io <https://sylabs.io/singularity/>`__ for more information on singularity images).
The only thing to note before you download this image is that the image is **~6.5Gb** so, depending on the internet speed, it will take **~15-30 min to download**.
The good news is that you should only need to do this once unless updates are made to the scripts or image.

Just download the singluarity image with:

  .. code-block:: bash

    wget https://www.dropbox.com/s/9g9k5u1xz23k4db/Demuxafy.sif
    wget https://www.dropbox.com/s/qx0q37aqa0j2gu9/Demuxafy.sif.md5



Then you should check to make sure that the image downloaded completely by comparing the image md5sum to the original md5sum.
You can do that by running the following commands:

  .. code-block:: bash

      md5sum Demuxafy.sif > downloaded_Demuxafy.sif.md5
      diff -s Demuxafy.sif.md5 downloaded_Demuxafy.sif.md5

If everything was downloaded correctly, that command should report:

  .. code-block:: bash

    Files Demuxafy.sif.md5 and downloaded_Demuxafy.sif.md5 are identical


.. note::

  Please note that the singularity image and this documentation is updated with each release. 
  This means that the most recent documentation may not be 100% compatible with the singularity image that you have.
  For example, additional parameters and functionality were implemented in v1.0.2 that was not available in v0.0.4.
  
  You can check the version of your singularity image to match with documentation with:

    .. code-block:: bash

      singularity inspect Demuxafy.sif


If you run into any issues with downloading the image or any issue with running anything from this image, you can reach out to us by submitting an issue at `Github <https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/issues>`__

.. admonition:: Demuxafy software versions - for the curious
  :class: dropdown

  Image build date: 24 July, 2022
 
    +----------------------------+---------------------------+-------------------------------+
    | Software Group             | Software                  | Version                       |
    +============================+===========================+===============================+
    |  Demultiplexing            | ``popscle``               |                               |
    |                            |  - ``demuxlet``           | v0.1-beta                     |
    |                            |  - ``freemuxlet``         |                               |
    |                            +---------------------------+-------------------------------+
    |                            | ``scSplit``               | v1.0.8.2                      |
    |                            +---------------------------+-------------------------------+
    |                            | ``Souporcell``            | v2.0                          |
    |                            +---------------------------+-------------------------------+
    |                            | ``Vireo``                 | v0.5.6                        |
    +----------------------------+---------------------------+-------------------------------+
    | Doublet Detecting          | ``DoubletDecon``          | v1.1.6                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``DoubletDetection``      | v3.0                          |
    |                            +---------------------------+-------------------------------+
    |                            | ``DoubletFinder``         | v2.0.3                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``scDblFinder``           | v1.8.0                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``scds``                  | v1.9.1                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``scrublet``              | v0.2.3                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``solo``                  | v1.2                          |
    +----------------------------+---------------------------+-------------------------------+
    | Supporting Softwares       | ``minimap2``              | v2.7-r654                     |
    |                            +---------------------------+-------------------------------+
    |                            | ``bedtools2``             | v2.28.0                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``vartrix``               | v1.1.3                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``htslib``                | v1.13 & v1.14                 |
    |                            +---------------------------+-------------------------------+
    |                            | ``samtools``              | v1.15                         |
    |                            +---------------------------+-------------------------------+
    |                            | ``bcftools``              | v1.13                         |
    |                            +---------------------------+-------------------------------+
    |                            | ``freebayes``             | v1.3.5                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``cellSNP-lite``          | v1.2.2                        |
    +----------------------------+---------------------------+-------------------------------+
    | R Supporting Packages      | ``argparse``              | v2.1.3                        |
    | (R v4.1.3)                 +---------------------------+-------------------------------+
    |                            | ``ComplexHeatmap``        | v2.10.0                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``ComplexUpset``          | v1.3.3                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``vcfR``                  | v1.12.0                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``Seurat``                | 4.1.0                         |
    |                            +---------------------------+-------------------------------+
    |                            | ``SingleCellExperiment``  | v1.16.0                       |
    +----------------------------+---------------------------+-------------------------------+
    | Python Supporting Packages | ``argparse``              | v1.4.0                        |
    | (Python v3.7.2)            +---------------------------+-------------------------------+
    |                            | ``numpy``                 | v1.21.2                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``matplotlib``            | v3.2.2                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``pandas``                | v1.3.5                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``PyVCF``                 | v0.6.8                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``scipy``                 | v1.7.3                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``scvi-tools``            | v0.14.6                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``umap-learn``            | v0.5.2                        |
    +----------------------------+---------------------------+-------------------------------+



              
