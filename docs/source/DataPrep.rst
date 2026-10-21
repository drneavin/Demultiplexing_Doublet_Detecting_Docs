.. _DataPrep-docs:

Data Preparation
================


There isn't a lot of data preparation to be done before running the demultiplexing or doublet detecting softwares.


Data Required
-------------
The demultiplexing and transcriptome-based doublet detecting softwares have different data input requirements:

+-------------------+-----------------------------------------------+------------------------------------------+
| Software Group    | .. centered:: Single Cell Count Data Required | .. centered:: SNP Genotype Data Required |
+===================+===============================================+==========================================+
| Demultiplexing    | .. centered:: |:heavy_check_mark:|            | .. centered:: |:heavy_check_mark:|       |
+-------------------+-----------------------------------------------+------------------------------------------+
| Doublet Detecting | .. centered:: |:heavy_check_mark:|            | .. centered:: |:heavy_multiplication_x:| |
+-------------------+-----------------------------------------------+------------------------------------------+


You won't need to pre-process the single cell count data unless you are using :ref:`DoubletFinder<doubletfinder-docs>` or :ref:`DoubletDecon<doubletdecon-docs>` which need QC-filtered and normalized counts (for example with `Seurat <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>`__).

For the demultiplexing softwares, you should filter the SNP genotypes that you will use.


SNP Genotype Data
---------------------------


.. NOTE::

  The SNP genotype data can be for multiplexed donors in the pool **OR** it can be publicly available common SNP genotypes.

  We provide instructions on how to access or prepare these data here.


Pulblicly Available SNP Genotype Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have provided common SNP vcf files that we have generated on both GRCh37 and GRCh38 and with a variety of filtering and 'chr' encoding.
These can be downloaded using the following links:

+------------------------+----------------------+--------------------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| .. centered:: Minor    | .. centered:: Genome | .. centered:: Region           | .. centered:: Chr Encoding   | .. centered:: vcf File                                                                                                                                            | .. centered:: md5sum File                                                                                                                                                 |
| .. centered:: Allele   |                      | .. centered:: Filtering        |                              |                                                                                                                                                                   |                                                                                                                                                                           |
| .. centered:: Frequency|                      |                                |                              |                                                                                                                                                                   |                                                                                                                                                                           |
+========================+======================+================================+==============================+===================================================================================================================================================================+===========================================================================================================================================================================+
| .. centered:: 1%       | .. centered:: GRCh37 | .. centered:: Genes            |  .. centered:: No 'chr'      | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf <https://drive.google.com/file/d/14_V1UciYvkAltR_Qwlo48YtXEeFaANGX/view?usp=sharing>`__                | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5 <https://drive.google.com/file/d/1WIFpnkfOnnixgwvO1U8FMjehnG1CfsfT/view?usp=sharing>`__                    |
|                        |                      |                                +------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      |                                | .. centered:: 'chr' encoding | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf <https://drive.google.com/file/d/12c2ktX7wu9QgFUAEuqWHBbX5bvb9e2e9/view?usp=sharing>`__          | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.md5 <https://drive.google.com/file/d/1FoIJeK7tgmbCxh8YZiWNdUyV1cHNmR0j/view?usp=sharing>`__              |
|                        |                      +--------------------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      | .. centered:: Exons            |  .. centered:: No 'chr'      | .. centered:: `GRCh37_1000G_MAF0.01_ExonFiltered_NoChr.vcf <https://drive.google.com/file/d/1-C---Bs1Xb_1_DF6H9rWtiz3_z1zI1R_/view?usp=sharing>`__                | .. centered:: `GRCh37_1000G_MAF0.01_ExonFiltered_NoChr.vcf.md5 <https://drive.google.com/file/d/1J7kahSWa2piWsNLSvRS21mmPIdECpKGz/view?usp=sharing>`__                    |
|                        |                      |                                +------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      |                                | .. centered:: 'chr' encoding | .. centered:: `GRCh37_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf <https://drive.google.com/file/d/1ojoR5Cn7Q3Pix1AdZVyIEgYtt3Nqs2bW/view?usp=drive_link>`__       | .. centered:: `GRCh37_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf.md5 <https://drive.google.com/file/d/1zWeqyLvvPf6O5OK0DfKuafraA_xyj67W/view?usp=sharing>`__              |
|                        +----------------------+--------------------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        | .. centered:: GRCh38 | .. centered:: Genes            |  .. centered:: No 'chr'      | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.vcf <https://drive.google.com/file/d/1i_yvsMVHwvmQAIrW7X3L283HrxQOgcZd/view?usp=sharing>`__                | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5 <https://drive.google.com/file/d/1h-iAwkplyfgEQAiuLnpBcBll-EfoIgby/view?usp=sharing>`__                    |
|                        |                      |                                +------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      |                                | .. centered:: 'chr' encoding | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf <https://drive.google.com/file/d/1WQ4qS85d0CHz1C2bS0qvPibx4x6duort/view?usp=sharing>`__          | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.md5 <https://drive.google.com/file/d/16L_-prwbv80wXGvviEodAJgd7CO_1Op-/view?usp=sharing>`__              |
|                        |                      +--------------------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      | .. centered:: Exons            |  .. centered:: No 'chr'      | .. centered:: `GRCh38_1000G_MAF0.01_ExonFiltered_NoChr.vcf <https://drive.google.com/file/d/184ZKLJDVqc9VSem0KWxj30dzsAfTANFI/view?usp=sharing>`__                | .. centered:: `GRCh38_1000G_MAF0.01_ExonFiltered_NoChr.vcf.md5 <https://drive.google.com/file/d/11mPnThwHcKeJ8qqGpvCaNl1q9QPVYri7/view?usp=sharing>`__                    |
|                        |                      |                                +------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      |                                | .. centered:: 'chr' encoding | .. centered:: `GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf <https://drive.google.com/file/d/1D8ZSHANlcmM9jVcvaBaymFBTt84yGr26/view?usp=sharing>`__          | .. centered:: `GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf.md5 <https://drive.google.com/file/d/111fRcd_lNpxA3pAERNT8H0-mwyAcUP8q/view?usp=sharing>`__              |
+------------------------+----------------------+--------------------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Of course, you can further filter these to a 5% minor allele frequency if you would prefer.

You can also download SNP genotype data and process it yourself from `1000G <https://www.internationalgenome.org/category/ftp/>`__ (hg19 and hg38) or `HRC <http://www.haplotype-reference-consortium.org/site>`__ (hg19 only).

For 1000G, use the instructions at the above link to access the data per your preferences and you can find the required files at the following directories:

  - The hg19 data is available at ``/ftp/release/``
  - The hg38 data is available at ``/ftp/release/20130502/supporting/GRCh38_positions/``


Preparing your own SNP Genotype Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is best to filter the SNP genotypes for common SNPs (generally > 1% or > 5% minor allele frequency) that overlap either exons or genes.
We typically suggest filtering for exons since it typically resultsin in ~250k SNPs to remain following filtering which is sufficient for demultiplexing without using too many SNPs which can slow down the demultiplexing softwares.
However, some capture types might be better suited to look for SNPs overlapping genes such as single nuclei RNA-seq.
For relative numbers of SNPs in the exons and introns, see the `issue raised by @jamesnemesh <https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/issues/49#issue-2182195018>`__.
Here we provide an example of how to do this filtering. 
We built the required softwares into the singularity image so you can run these filtering steps with the image.

.. NOTE::
  
  We have found it best to impute reference SNP genotypes so there are more SNP locations available. 
  If you are using reference SNP genotypes for the donors in your pool, please be sure to impute before filtering.

Filter for Common SNPs
^^^^^^^^^^^^^^^^^^^^^^
First, filter the SNP genotypes for common SNPs - 5% minor allele frequency should work for most datasets but you can change this to another minor allele frequency if you would like.
We assume that you have already filtered imputed SNPs for quality based on the imputation method used.

.. code-block:: bash

  singularity exec Demuxafy.sif bcftools filter --include 'MAF>=0.05' -Oz --output $OUTDIR/common_maf0.05.vcf.gz $VCF

Where ``$OUTDIR`` is the output directory where you want to save the results and ``$VCF`` is the path to the SNP genotype vcf file.

Filter for SNPs overlapping Exons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Next, filter for the SNPs that overlap exons.
The below command uses `vcftools`` with a `bed file <https://en.wikipedia.org/wiki/BED_(file_format)>` that contains the locations of each exon (or gene if you prefer) across the genome. 


.. NOTE::

  You can get an exon bed using the `UCSC table browser <https://genome.ucsc.edu/cgi-bin/hgTables>`__ (see instructions `here <https://www.biostars.org/p/93011/>`__) and we have also provided bed files for :download:`hg19 <../../references/hg19exonsUCSC.bed>` and :download:`hg38 <../../references/hg38exonsUCSC.bed>`

  Here is an example of what the head of a bed file might look like and you can find a description of the file type `here <https://en.wikipedia.org/wiki/BED_(file_format)>`.
  The only required columns are the first three that contain the location of the exon:

  .. code-block: bash

    1       11868   12227   ENST00000456328.2_exon_0_0_chr1_11869_f 0       +
    1       12612   12721   ENST00000456328.2_exon_1_0_chr1_12613_f 0       +
    1       13220   14409   ENST00000456328.2_exon_2_0_chr1_13221_f 0       +
    1       12009   12057   ENST00000450305.2_exon_0_0_chr1_12010_f 0       +
    1       12178   12227   ENST00000450305.2_exon_1_0_chr1_12179_f 0       +
    1       12612   12697   ENST00000450305.2_exon_2_0_chr1_12613_f 0       +
    1       12974   13052   ENST00000450305.2_exon_3_0_chr1_12975_f 0       +
    1       13220   13374   ENST00000450305.2_exon_4_0_chr1_13221_f 0       +
    1       13452   13670   ENST00000450305.2_exon_5_0_chr1_13453_f 0       +
    1       14403   14501   ENST00000488147.1_exon_0_0_chr1_14404_r 0       -

.. code-block:: bash

  singularity exec Demuxafy.sif vcftools \
    --gzvcf $OUTDIR/common_maf0.05.vcf.gz \
    --max-alleles 2 \
    --remove-indels \
    --bed $BED \
    --recode \
    --recode-INFO-all \
    --out $OUTDIR/common_maf0.05_exon_filtered


We typically expect ~250k SNPs to remain following this filtering step.


.. _TestData:

Test Dataset
------------
In addition, we have provided a test dataset that you can use.
Find the data that can be downloaded below along with some information about the data.
The data have been aligned to GRCh38 so for running.
If you are running souporcell, you will needa a fasta file, which you can download from `Ensembl FTP <https://asia.ensembl.org/info/data/ftp/index.html>`__

.. admonition:: Information
  :class: important

    The test dataset includes 20,982 droplets captured of PBMCs from 13 multiplexed individuals.

10x Directories + Other Necessary Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We have provided this dataset as the complete dataset which is pretty large (~40Gb tar.gz directory).
Therefore, we have also provided the same dataset where the data has been significantly reduced.

.. WARNING:: 
	The reduced test dataset may not produce real-world results due to the small size - especially for doublet detecting softwares since the reads have been significantly downsampled to reduce the size.

You can download the dataset with one of the following commands:

.. tabs::

  .. tab:: Complete Dataset

    First, download the dataset and the md5sum:

    .. code-block:: bash

      wget -O TestData4PipelineFull.tar.gz 'https://drive.usercontent.google.com/download?id=1RzEYVEvhpWvjwjs671UyTV3R8iO4VbwC&export=download&confirm=t'
      wget -O TestData4PipelineFull.tar.gz.md5 'https://drive.usercontent.google.com/download?id=https://drive.google.com/file/d/1umG5tOSJwpwy9iIVp7jCjoFT7wDU9ctZ/view?usp=sharing&export=download&confirm=t'

    After downloading the tar.gz directory, it is best to make sure the md5sum of the ``TestData4PipelineFull.tar.gz`` file matches the md5sum in the ``TestData4PipelineFull.tar.gz.md5``:

    .. code-block:: bash

      md5sum TestData4PipelineFull.tar.gz > downloaded_TestData4PipelineFull.tar.gz.md5
      diff -s TestData4PipelineFull.tar.gz.md5 downloaded_TestData4PipelineFull.tar.gz.md5

    That should return the following statement indicating that the two md5sums are identical:

    .. code-block:: bash

      Files TestData4PipelineFull.tar.gz.md5 and downloaded_TestData4PipelineFull.tar.gz.md5 are identical

    Finally, you can access the data by unzipping the file:

    .. code-block:: bash

      tar -xvf TestData4PipelineFull.tar.gz

    This should unzip the ``TestData4PipelineFull`` directory where you will have the following file structure:

    .. code-block:: bash

      TestData4PipelineFull
      ├── donor_list.txt 
      ├── samplesheet.txt
      ├── test_dataset
      │   ├── outs
      │   │   └── filtered_gene_bc_matrices
      │   │       └── Homo_sapiens_GRCh38p10
      │   │           ├── barcodes.tsv
      │   │           ├── genes.tsv
      │   │           └── matrix.mtx
      │   ├── possorted_genome_bam.bam
      │   └── possorted_genome_bam.bam.bai
      └── test_dataset.vcf


  .. tab:: Reduced Dataset

    First, download the reduced dataset and the md5sum:

    .. code-block:: bash

      wget -O TestData4PipelineSmall.tar.gz 'https://drive.usercontent.google.com/download?id=1cx2Wv46hoYwn659ixSFcYGVKKQw8EXDy&export=download&confirm=t'
      wget -O TestData4PipelineSmall.tar.gz.md5 'https://drive.usercontent.google.com/download?id=13zuR242x49xGbiZikbnaySA8DgzA1P-W&export=download&confirm=t'

    After downloading the tar.gz directory, it is best to make sure the md5sum of the ``TestData4PipelineSmall.tar.gz`` file matches the md5sum in the ``TestData4PipelineSmall.tar.gz.md5``:

    .. code-block:: bash

      md5sum TestData4PipelineSmall.tar.gz > downloaded_TestData4PipelineSmall.tar.gz.md5
      diff -s TestData4PipelineSmall.tar.gz.md5 downloaded_TestData4PipelineSmall.tar.gz.md5

    That should return the following statement indicating that the two md5sums are identical:

    .. code-block:: bash

      Files TestData4PipelineSmall.tar.gz.md5 and downloaded_TestData4PipelineSmall.tar.gz.md5 are identical

    Finally, you can access the data by unzipping the file:

    .. code-block:: bash

      tar -xvf TestData4PipelineSmall.tar.gz

    This should unzip the ``TestData4PipelineSmall`` directory where you will have the following file structure:

    .. code-block:: bash

      TestData4PipelineSmall
      ├── donor_list.txt
      ├── individuals_list_dir
      │   └── test_dataset.txt
      ├── samplesheet.txt
      ├── test_dataset
      │   └── outs
      │       ├── filtered_gene_bc_matrices
      │       │   └── Homo_sapiens_GRCh38p10
      │       │       ├── barcodes.tsv
      │       │       ├── genes.tsv
      │       │       └── matrix.mtx
      │       ├── pooled.sorted.bam
      │       └── pooled.sorted.bam.bai
      └── test_dataset.vcf

Seurat Object
^^^^^^^^^^^^^^
We have also provided a filtered, QC normalized Seurat object (needed for :ref:`DoubletFinder<doubletfinder-docs>` and :ref:`DoubletDecon<doubletdecon-docs>`)

Download the rds object and the md5sum:

.. code-block:: bash

  wget -O TestData_Seurat.rds 'https://drive.usercontent.google.com/download?id=1OL5jDyfvzBkCYn7V_GMGj-wjhYqmnYpx&export=download&confirm=t'
  wget -O TestData_Seurat.rds.md5 'https://drive.usercontent.google.com/download?id=1BRFJrs03zZPqSBwYq5oM11TUkRolAKFR&export=download&confirm=t'


After downloading the rds.object, it is best to make sure the md5sum of the ``TestData_Seurat.rds`` file matches the md5sum in the ``TestData_Seurat.rds.md5``:

.. code-block:: bash

	md5sum TestData_Seurat.rds > downloaded_TestData_Seurat.rds.md5
	diff -s TestData_Seurat.rds.md5 downloaded_TestData_Seurat.rds.md5

That should return the following statement indicating that the two md5sums are identical:

.. code-block:: bash

	Files TestData_Seurat.rds.md5 and downloaded_TestData_Seurat.rds.md5 are identical

The ``TestData_Seurat.rds`` can then be used directly as input for the :ref:`DoubletFinder<doubletfinder-docs>` and :ref:`DoubletDecon<doubletdecon-docs>` tutorials.
You can also load the ``TestData_Seurat.rds`` into ``R`` to see the seurat object by first opening R:

.. code-block:: bash

  singularity exec Demuxafy.sif R

Then read it in with:

.. code-block:: R

  seurat <- readRDS(TestData_Seurat.rds)


.. Note:: 
  We have used this dataset for each of the tutorials.
  The example tables in the *Results and Interpretation* sections of each tutorial are the results from this dataset.


