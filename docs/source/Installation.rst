Installation
==========================

Singularity Image Download
--------------------------------
Installation should be pretty painless (we hope).
We have  provided all the softwares in a singularity image which provides continuity across different computing platforms (see `HPCNG Singluarity <https://singularity.hpcng.org/>`__ and `Sylabs io <https://sylabs.io/singularity/>`__ for more information on singularity images).
The only thing to note before you download this image is that the image is **~7.5Gb** so, depending on the internet speed, it will take **~15-30 min to download**.
The good news is that you should only need to do this once unless updates are made to the scripts or image.

Just download the singlularity image with:

  .. code-block:: bash

    wget -O Demuxafy.sif 'https://www.dropbox.com/scl/fi/kykwi78vk4yifbbag5ajz/Demuxafy.sif?rlkey=5hcugu6ztpy0eik3xno63xiar'
    wget -O Demuxafy.sif.md5 'https://www.dropbox.com/scl/fi/37oj9y1frzhqazl4h8s21/Demuxafy.sif.md5?rlkey=o2bn5wp9q68numlaav8gg95kh'




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
  There is documentation to match each release so be sure you're using documentation that matches your singularity image.
  
  You can check the version of your singularity image to match with documentation with:

    .. code-block:: bash

      singularity inspect Demuxafy.sif


If you run into any issues with downloading the image or any issue with running anything from this image, you can reach out to us by submitting an issue at `Github <https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/issues>`__

.. admonition:: Demuxafy software versions - for the curious
  :class: dropdown

  Image version: 3.0.0
  Image build date: 5 June, 2024
 
    +----------------------------+---------------------------+-------------------------------+
    | Software Group             | Software                  | Version                       |
    +============================+===========================+===============================+
    | Demultiplexing             | ``Demuxalot``             | v0.4.0                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``popscle``               |                               |
    |                            |  - ``demuxlet``           | v0.1-beta                     |
    |                            |  - ``freemuxlet``         |                               |
    |                            +---------------------------+-------------------------------+
    |                            | ``Dropulation``           | v3.0.1-7-gf7acbbc             |
    |                            +---------------------------+-------------------------------+
    |                            | ``scSplit``               | v1.0.9                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``Souporcell``            | v2.5                          |
    |                            +---------------------------+-------------------------------+
    |                            | ``Vireo``                 | v0.5.8                        |
    +----------------------------+---------------------------+-------------------------------+
    | Doublet Detecting          | ``DoubletDecon``          | v1.1.6                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``DoubletDetection``      | v4.2                          |
    |                            +---------------------------+-------------------------------+
    |                            | ``DoubletFinder``         | v2.0.4                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``scDblFinder``           | v1.18.0                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``scds``                  | v1.13.1                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``scrublet``              | v0.2.3                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``solo``                  | v1.2                          |
    +----------------------------+---------------------------+-------------------------------+
    | Supporting Softwares       | ``minimap2``              | 2.26-r1175                    |
    |                            +---------------------------+-------------------------------+
    |                            | ``bedtools2``             | v2.30.0                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``vartrix``               | v1.1.22                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``htslib``                | v1.16                         |
    |                            +---------------------------+-------------------------------+
    |                            | ``samtools``              | v1.20                         |
    |                            +---------------------------+-------------------------------+
    |                            | ``bcftools``              | v1.20                         |
    |                            +---------------------------+-------------------------------+
    |                            | ``freebayes``             | v1.3.7                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``cellSNP-lite``          | v1.2.3                        |
    +----------------------------+---------------------------+-------------------------------+
    | R Supporting Packages      | ``argparse``              | v2.2.3                        |
    | (R v4.4.0)                 +---------------------------+-------------------------------+
    |                            | ``ComplexHeatmap``        | v2.20.0                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``ComplexUpset``          | v1.3.3                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``vcfR``                  | v1.15.0                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``Seurat``                | v5.1.0                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``SingleCellExperiment``  | v1.26.0                       |
    +----------------------------+---------------------------+-------------------------------+
    | Python Supporting Packages | ``argparse``              | v1.1.0                        |
    | (Python v3.9.7)            +---------------------------+-------------------------------+
    |                            | ``numpy``                 | v1.25.2                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``matplotlib``            | v3.9.0                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``pandas``                | v1.5.3                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``PyVCF``                 | v0.6.8                        |
    |                            +---------------------------+-------------------------------+
    |                            | ``scipy``                 | v1.10.1                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``scvi-tools``            | v0.14.6                       |
    |                            +---------------------------+-------------------------------+
    |                            | ``umap-learn``            | v0.5.6                        |
    +----------------------------+---------------------------+-------------------------------+



              

.. _Singularity-docs:


Notes About Singularity Images
--------------------------------

Singularity images effectively store an operating system with files, softwares etc. that can be easily transported across different operating systems - ensuring reproducibility.
Most HPCs have singularity installed making it easy to implement.
There are some tips and tricks we have identified through using singularity images that we thought might help new users.

Tips and Tricks
++++++++++++++++++
1. Error: File Not Found
^^^^^^^^^^^^^^^^^^^^^^^^
  **Reason**

  Singularity only loads the directories directly downstream from where you execute the singularity command.
  If any of the files that need to be accessed by the command are not downstream of the that location, you will receive an error similar to this one:

  .. code-block:: bash

    Failed to open file "/path/to/readfile.tsv" : No such file or directory

  If you then check for that file:

  .. code-block:: bash

    ll /path/to/readfile.tsv

  We can see that the  file does truly exist:

  .. code-block:: bash

    -rw-rw-r-- 1 user group 70636291 Dec 21  2020 /path/to/readfile.tsv

  **Solution**

  The easiest solution to this problem is to "bind" a path upstream of all the files that will need to be accessed by your command:

  .. code-block:: bash

    singularity exec --bind /path Demuxafy.sif ...


If you don't have access to Singularity on your HPC, you can ask your HPC administrators to install it (see the `Singularity page <https://sylabs.io/guides/3.0/user-guide/quick_start.html>`__)
