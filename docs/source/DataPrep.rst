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

.. NOTE::

  The SNP genotype data can be for multiplexed donors in the pool **OR** it can be publicly available common SNP genotypes which can be downloaded from 1000G (`hg19 <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/>`__ or `hg38 <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/>`__ or from HRC (`hg19 <http://www.haplotype-reference-consortium.org/site>`__.


There isn't any pre-processing that's needed for the single cell count data.
For the demultiplexing softwares, you should filter the SNP genotypes that you will use.
It is best to filter the SNP genotypes for common SNPs (generally > 1% or > 5% minor allele frequency) that overlap exons.
Here we provide an example of how to do this filtering. 
We built the required softwares into the singularity image so you can run the filtering with the image.

.. NOTE::
  
  We have found it best to impute reference SNP genotypes so there are more SNP locations available. 
  If you are using reference SNP genotypes for the donors in your pool, please be sure to impute before filtering.

First, filter the imputed reference SNP genotypes or the publicly available SNP genotypes for common SNPs - here we use 5% minor allele frequency.

.. code-block:: bash

  singularity exec Demuxafy.sif bcftools filter --include 'MAF>=0.05' -Oz --output $OUTDIR/common_maf0.05.vcf.gz $VCF

Next, filter for the SNPs that overlap exons.

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
In addition, we have provided test data that you can use to ensure the singularity image is working correctly.

.. admonition:: Information
  :class: important

    The test dataset includes 20,982 droplets captured of PBMCs from 14 multiplexed individuals.

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

.. Note:: 
  We have used this dataset for each of the tutorials.
  The example tables in the *Results and Interpretation* sections of each tutorial are the results from this dataset.
