.. _Dropulation-docs:


Dropulation
===========================

.. _Dropulation: https://github.com/broadinstitute/Drop-seq/blob/master/doc/Census-seq_Computational_Protcools.pdf
.. _preprint: https://www.biorxiv.org/content/10.1101/2022.03.07.483367v1
.. _GENCODE: https://www.gencodegenes.org/human/



Dropulation_ is a genotype demultiplexing software that requires reference genotypes to be available for each individual in the pool. 
Therefore, if you don't have reference genotypes, you may want to demultiplex with one of the softwares that do not require reference genotype data
(:ref:`Freemuxlet <Freemuxlet-docs>`, :ref:`scSplit <scSplit-docs>`, :ref:`Souporcell <Souporcell-docs>` or :ref:`Vireo <Vireo-docs>`)




Data
----
This is the data that you will need to have prepare to run Dropulation_:

.. admonition:: Required
  :class: important

  - Reference SNP genotypes for each individual (``$VCF``)

    - Filter for common SNPs (> 5% minor allele frequency) and SNPs overlapping genes

  - Barcode file (``$BARCODES``)

  - Bam file (``$BAM``)

    - Aligned single cell reads

    - This must contain gene annotations as produced by ``TagReadWithGeneFunction`` to have gene name (gn), gene strand (gs) and gene function (gf) - we have included code for doing this below, it will require a gtf file for annotation.

  - Output directory (``$DROPULATION_OUTDIR``)

  - GTF gene annotation file (``$GTF``)

    - If you already have a bam annotated with ``TagReadWithGeneFunction``, then you won't need the ``$GTF`` file

    - Ideally, this would be the gtf file used for alignment but you can also get gtf files from GENCODE_

  - A text file with the individual ids (``$INDS``)
  
    - File containing the individual ids (separated by line) as they appear in the vcf file

    - For example, this is the :download:`individual file <_download_files/Individuals.txt>` for our example dataset


Run Dropulation
-----------------
First, let's assign the variables that will be used to execute each step.

.. admonition:: Example Variable Settings
  :class: grey

    Below is an example of the variables that we can set up to be used in the command below.
    These are files provided as a :ref:`test dataset <TestData>` available in the :ref:`Data Preparation Documentation <DataPrep-docs>`
    Please replace paths with the full path to your data on your system.

    .. code-block:: bash

      VCF=/path/to/TestData4PipelineFull/test_dataset.vcf
      BARCODES=/path/to/TestData4PipelineFull/test_dataset/outs/filtered_gene_bc_matrices/Homo_sapiens_GRCh38p10/barcodes.tsv
      BAM=/path/to/test_dataset/possorted_genome_bam.bam
      DROPULATION_OUTDIR=/path/to/output/dropulation
      INDS=/path/to/TestData4PipelineFull/donor_list.txt
      GTF=/path/to/genes.gtf


Bam Annotation
^^^^^^^^^^^^^^^^^^^^
You will most likely need to annotate your bam using ``TagReadWithGeneFunction`` (unless you have already annotated your bam with this function).
Please note that the ``\`` at the end of each line is purely for readability to put a separate parameter argument on each line.

  .. code-block:: bash

    singularity exec Demuxafy.sif TagReadWithGeneFunction \
              --ANNOTATIONS_FILE $GTF \
              --INPUT $BAM \
              --OUTPUT $DROPULATION_OUTDIR/possorted_genome_bam_dropulation_tag.bam


If the bam annotation is successful, you will have these new files in your ``$DROPULATION_OUTDIR``:

.. code-block:: bash

  /path/to/output/dropulation
  └── possorted_genome_bam_dropulation_tag.bam



Dropulation Assignment
^^^^^^^^^^^^^^^^^^^^^^^^^
First, we will identify the most likely singlet donor for each droplet.

.. admonition:: Note
  :class: note

  Please change the cell barcode and molecular barcode tags as necessary. 
  For 10x experiments processed with cellranger, this should be 'CB' for the ``CELL_BARCODE_TAG`` and 'UB' for the ``MOLECULAR_BARCODE_TAG``

Please note that the ``\`` at the end of each line is purely for readability to put a separate parameter argument on each line.

.. code-block:: bash

  AssignCellsToSamples --CELL_BC_FILE $BARCODES \
            --INPUT_BAM $DROPULATION_OUTDIR/possorted_genome_bam_dropulation_tag.bam \
            --OUTPUT $DROPULATION_OUTDIR/assignments.tsv.gz \
            --VCF $VCF \
            --SAMPLE_FILE $INDS \
            --CELL_BARCODE_TAG 'CB' \
            --MOLECULAR_BARCODE_TAG 'UB' \
            --VCF_OUTPUT $DROPULATION_OUTDIR/assignment.vcf \
            --MAX_ERROR_RATE 0.05


If the bam annotation is successful, you will have these new files in your ``$DROPULATION_OUTDIR``:

.. code-block:: bash
  :emphasize-lines: 2,3,4

  /path/to/output/dropulation
  ├── assignments.tsv.gz
  ├── out_vcf.vcf
  ├── out_vcf.vcf.idx
  └── possorted_genome_bam_dropulation_tag.bam


Dropulation Doublet
^^^^^^^^^^^^^^^^^^^^^^^^^
Next, we will identify the likelihoods of each droplet being a doublet.

.. admonition:: Note
  :class: note

  Please change the cell barcode and molecular barcode tags as necessary. 
  For 10x experiments processed with cellranger, this should be 'CB' for the ``CELL_BARCODE_TAG`` and 'UB' for the ``MOLECULAR_BARCODE_TAG``

Please note that the ``\`` at the end of each line is purely for readability to put a separate parameter argument on each line.

.. code-block:: bash

  DetectDoublets --CELL_BC_FILE $BARCODES \
            --INPUT_BAM $DROPULATION_OUTDIR/possorted_genome_bam_dropulation_tag.bam \
            --OUTPUT $DROPULATION_OUTDIR/likelihoods.tsv.gz \
            --VCF $VCF \
            --CELL_BARCODE_TAG 'CB' \
            --MOLECULAR_BARCODE_TAG 'UB' \
            --SINGLE_DONOR_LIKELIHOOD_FILE $DROPULATION_OUTDIR/assignments.tsv.gz \
            --SAMPLE_FILE $INDS \
            --MAX_ERROR_RATE 0.05


Dropulation Call
^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, we will make final assignments for each droplet based on the doublet and assignment calls.

Please note that the ``\`` at the end of each line is purely for readability to put a separate parameter argument on each line.

.. code-block:: bash

  Rscript dropulation_call.R --assign $DROPULATION_OUTDIR/assignments.tsv.gz \
                             --doublet $DROPULATION_OUTDIR/likelihoods.tsv.gz \
                             --out $DROPULATION_OUTDIR/updated_assignments.tsv.gz


If the bam annotation is successful, you will have these new files in your ``$DROPULATION_OUTDIR``:

.. code-block:: bash
  :emphasize-lines: 5

  /path/to/output/dropulation
  ├── assignments.tsv.gz
  ├── out_vcf.vcf
  ├── out_vcf.vcf.idx
  ├── possorted_genome_bam_dropulation_tag.bam
  └── updated_assignments.tsv.gz
  


Dropulation Summary
^^^^^^^^^^^^^^^^^^^^^
We have provided a script that will summarize the number of droplets classified as doublets, ambiguous and assigned to each donor by Dropulation_ and write it to the ``$DROPULATION_OUTDIR``. 
You can run this to get a fast and easy summary of your results by providing the path to your result file:

.. code-block:: bash

  singularity exec Demuxafy.sif bash Dropulation_summary.sh $DROPULATION_OUTDIR/updated_assignments.tsv.gz


which will return:

  +-----------------+--------------+
  | Classification  | Assignment N |
  +=================+==============+
  | 113_113         | 1334         |
  +-----------------+--------------+
  | 349_350         | 1458         |
  +-----------------+--------------+
  | 352_353         | 1607         |
  +-----------------+--------------+
  | 39_39           | 1297         |
  +-----------------+--------------+
  | 40_40           | 1078         |
  +-----------------+--------------+
  | 41_41           | 1127         |
  +-----------------+--------------+
  | 42_42           | 1419         |
  +-----------------+--------------+
  | 43_43           | 1553         |
  +-----------------+--------------+
  | 465_466         | 1094         |
  +-----------------+--------------+
  | 596_597         | 1255         |
  +-----------------+--------------+
  | 597_598         | 1517         |
  +-----------------+--------------+
  | 632_633         | 868          |
  +-----------------+--------------+
  | 633_634         | 960          |
  +-----------------+--------------+
  | 660_661         | 1362         |
  +-----------------+--------------+
  | doublet         | 3053         |
  +-----------------+--------------+



or you can write it straight to a file:

.. code-block:: bash

  singularity exec Demuxafy.sif bash Dropulation_summary.sh $DROPULATION_OUTDIR/updated_assignments.tsv.gz > $DROPULATION_OUTDIR/dropulation_summary.tsv


.. admonition:: Note

  To check if these numbers are consistent with the expected doublet rate in your dataset, you can use our `Doublet Estimation Calculator <test.html>`__.



.. _dropulation-results:

Dropulation Results and Interpretation
----------------------------------------
After running the Dropulation_ steps and summarizing the results, you will have a number of files from some of the intermediary steps. 
These are the files that most users will find the most informative:

  - ``updated_assignments.tsv.gz``

    - The predicted annotations for each droplet and metrics:
      
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      | Barcode                    | dropulation_Likelihood | dropulation_Assignment | dropulation_DropletType | dropulation_Nsnps |      dropulation_Numis  |
      +============================+========================+========================+=========================+===================+=========================+
      |                            |                        |                        |                         |                   |                         |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      |                            |                        |                        |                         |                   |                         |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      |                            |                        |                        |                         |                   |                         |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      |                            |                        |                        |                         |                   |                         |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      |                            |                        |                        |                         |                   |                         |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      | ...                        | ...                    | ...                    |  ...                    | ...               | ...                     |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+



Merging Results with Other Software Results
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.

Citation
-----------
If you used the Demuxafy platform for analysis, please reference our preprint_ as well as Dropulation_.




