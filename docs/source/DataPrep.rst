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
| .. centered:: 1%       | .. centered:: GRCh37 | .. centered:: Genes            |  .. centered:: No 'chr'      | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf <https://www.dropbox.com/s/rce452zoawd0eee/GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf>`__             | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5 <https://www.dropbox.com/s/jrfb287hux6ehtg/GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5>`__             |
|                        |                      |                                +------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      |                                | .. centered:: 'chr' encoding | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf <https://www.dropbox.com/s/qrn6df1i1wxukxn/GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf>`__ | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.md5 <https://www.dropbox.com/s/on47saot3d2cgij/GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.md5>`__ |
|                        |                      +--------------------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      | .. centered:: Exons            |  .. centered:: No 'chr'      | .. centered:: `GRCh37_1000G_MAF0.01_ExonFiltered_NoChr.vcf <https://www.dropbox.com/s/fcempba0u0xwz2q/GRCh37_1000G_MAF0.01_ExonFiltered_NoChr.vcf>`__             | .. centered:: `GRCh37_1000G_MAF0.01_ExonFiltered_NoChr.vcf.md5 <https://www.dropbox.com/s/t4kucrvg4h3bkpj/GRCh37_1000G_MAF0.01_ExonFiltered_NoChr.vcf.md5>`__             |
|                        |                      |                                +------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      |                                | .. centered:: 'chr' encoding | .. centered:: `GRCh37_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf <https://www.dropbox.com/s/qa6qbxoodwrpl32/GRCh37_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf>`__ | .. centered:: `GRCh37_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf.md5 <https://www.dropbox.com/s/0ra5c19mta40aq9/GRCh37_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf.md5>`__ |
|                        +----------------------+--------------------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        | .. centered:: GRCh38 | .. centered:: Genes            |  .. centered:: No 'chr'      | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.vcf <https://www.dropbox.com/s/4nmm344g4j7pou4/GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.vcf>`__             | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5 <https://www.dropbox.com/s/izwp3l8oqwrt9dn/GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5>`__             |
|                        |                      |                                +------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      |                                | .. centered:: 'chr' encoding | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf <https://www.dropbox.com/s/ycfxs407sqgoori/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf>`__ | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.md5 <https://www.dropbox.com/s/77opxja1cgaq994/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.md5>`__ |
|                        |                      +--------------------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      | .. centered:: Exons            |  .. centered:: No 'chr'      | .. centered:: `GRCh38_1000G_MAF0.01_ExonFiltered_NoChr.vcf <https://www.dropbox.com/s/myyhcvhlxb0touq/GRCh38_1000G_MAF0.01_ExonFiltered_NoChr.vcf>`__             | .. centered:: `GRCh38_1000G_MAF0.01_ExonFiltered_NoChr.vcf.md5 <https://www.dropbox.com/s/yv9johlb3ogrbks/GRCh38_1000G_MAF0.01_ExonFiltered_NoChr.vcf.md5>`__             |
|                        |                      |                                +------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                        |                      |                                | .. centered:: 'chr' encoding | .. centered:: `GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf <https://www.dropbox.com/s/1akkxa7s9u485ka/GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf>`__ | .. centered:: `GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf.md5 <https://www.dropbox.com/s/aobteqk3z820ckn/GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding.vcf.md5>`__ |
+------------------------+----------------------+--------------------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Of course, you can further filter these to a 5% minor allele frequency if you would prefer.

You can also download SNP genotype data and process it yourself from `1000G <https://www.internationalgenome.org/category/ftp/>`__ (hg19 and hg38) or `HRC <http://www.haplotype-reference-consortium.org/site>`__ (hg19 only).

For 1000G, use the instructions at the above link to access the data per your preferences and you can find the required files at the following directories:

  - The hg19 data is available at ``/ftp/release/``
  - The hg38 data is available at ``/ftp/release/20130502/supporting/GRCh38_positions/``


Preparing your own SNP Genotype Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is best to filter the SNP genotypes for common SNPs (generally > 1% or > 5% minor allele frequency) that overlap exons.
Here we provide an example of how to do this filtering. 
We built the required softwares into the singularity image so you can run these filtering steps with the image.

.. NOTE::
  
  We have found it best to impute reference SNP genotypes so there are more SNP locations available. 
  If you are using reference SNP genotypes for the donors in your pool, please be sure to impute before filtering.

Filter for Common SNPs
^^^^^^^^^^^^^^^^^^^^^^
First, filter the SNP genotypes for common SNPs - 5% minor allele frequency should work for most datasets but you can change this to another minor allele frequency if you would like.

.. code-block:: bash

  singularity exec Demuxafy.sif bcftools filter --include 'MAF>=0.05' -Oz --output $OUTDIR/common_maf0.05.vcf.gz $VCF

Where ``$OUTDIR`` is the output directory where you want to save the results and ``$VCF`` is the path to the SNP genotype vcf file.

Filter for SNPs overlapping Exons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Next, filter for the SNPs that overlap exons.

.. NOTE::

  You can get an exon bed using the `UCSC table browser <https://genome.ucsc.edu/cgi-bin/hgTables>`__ (see instructions `here <https://www.biostars.org/p/93011/>`__) and we have also provided bed files for :download:`hg19 <../../references/hg19exonsUCSC.bed>` and :download:`hg38 <../../references/hg38exonsUCSC.bed>`

.. code-block:: bash

  singularity exec Demuxafy.sif vcftools \
    --gzvcf $OUTDIR/common_maf0.05.vcf.gz \
    --max-alleles 2 \
    --remove-indels \
    --bed $BED \
    --recode \
    --recode-INFO-all \
    --out $OUTDIR/common_maf0.05_exon_filtered



Test Dataset
------------
In addition, we have provided a test dataset that you can use.

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

      wget https://www.dropbox.com/s/3oujqq98y400rzz/TestData4PipelineFull.tar.gz
      wget https://www.dropbox.com/s/5n7u723okkf5m3l/TestData4PipelineFull.tar.gz.md5

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

    This should unzip the `TestData4PipelineFull` directory where you will have the following file structure:

    .. code-block:: bash

      TestData4PipelineFull
      ├── donor_list.txt
      ├── individuals_list_dir
      │   └── test_dataset.txt
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

      wget https://www.dropbox.com/s/m8u61jn4i1mcktp/TestData4PipelineSmall.tar.gz
      wget https://www.dropbox.com/s/ykjg86q3xw39wqr/TestData4PipelineSmall.tar.gz.md5

    After downloading the tar.gz directory, it is best to make sure the md5sum of the ``TestData4PipelineSmall.tar.gz`` file matches the md5sum in the ``TestData4PipelineSmall.tar.gz.md5``:

    .. code-block:: bash

      md5sum TestData4PipelineSmall.tar.gz > downloaded_TestData4PipelineSmall.tar.gz.md5
      diff -s TestData4PipelineSmall.tar.gz.md5 downloaded_TestData4PipelineSmall.tar.gz.md5

    That should return the following statement indicating that the two md5sums are identical:

    .. code-block:: bash

      Files TestData4PipelineSmall.tar.gz.md5 and downloaded_TestData4PipelineSmall.tar.gz.md5 are identical

    Finally, you can access the data by unzipping the file:

    .. code-block:: bash

      tar -xvf TestData4PipelineFull.tar.gz

    This should unzip the `TestData4PipelineSmall` directory where you will have the following file structure:

    .. code-block:: bash

      TestData4PipelineSmall
      ├── donor_list.txt
      ├── individuals_list_dir
      │   └── test_dataset.txt
      ├── ls
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

	wget https://www.dropbox.com/s/po4gy2j3eqohhjv/TestData_Seurat.rds
	wget https://www.dropbox.com/s/rmix7tt9aw28n7i/TestData_Seurat.rds.md5


After downloading the rds.object, it is best to make sure the md5sum of the ``TestData_Seurat.rds`` file matches the md5sum in the ``TestData_Seurat.rds.md5``:

.. code-block:: bash

	md5sum TestData_Seurat.rds > downloaded_TestData_Seurat.rds.md5
	diff -s TestData_Seurat.rds.md5 downloaded_TestData_Seurat.rds.md5

That should return the following statement indicating that the two md5sums are identical:

.. code-block:: bash

	Files TestData_Seurat.rds.md5 and downloaded_TestData_Seurat.rds.md5 are identical

The `TestData_Seurat.rds` can then be used directly as input for the :ref:`DoubletFinder<doubletfinder-docs>` and :ref:`DoubletDecon<doubletdecon-docs>` tutorials.
You can also load the `TestData_Seurat.rds` into `R` to see the seurat object by first opening R:

.. code-block:: bash

  singularity exec Demuxafy.sif R

Then read it in with:

.. code-block:: R

  seurat <- readRDS(TestData_Seurat.rds)


.. Note:: 
  We have used this dataset for each of the tutorials.
  The example tables in the *Results and Interpretation* sections of each tutorial are the results from this dataset.


