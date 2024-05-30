.. _Dropulation-docs:


Dropulation
===========================

.. _Dropulation: https://github.com/broadinstitute/Drop-seq/blob/master/doc/Census-seq_Computational_Protcools.pdf
.. _publication: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03224-8
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


.. admonition:: Optional

    - The SAM tag used in the Bam file to annotate the aligned single cell reads with their corresponding cell barcode (``$CELL_TAG``)

      - If not specified, _Demuxlet defaults to using ``XC``.

    - The SAM tag used in the Bam file to annotate the aligned single cell reads with their corresponding unique molecular identifier (UMI) (``$UMI_TAG``)

      - If not specified, _Demuxlet defaults to using ``XM``.


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
      CELL_TAG=CB  # default: XC
      UMI_TAG=UB  # default: XM


Bam Annotation
^^^^^^^^^^^^^^^^^^^^
.. admonition:: |:stopwatch:| Expected Resource Usage
  :class: note

  ~4h using a total of 3Gb memory when using 12 thread for the full :ref:`Test Dataset <TestData>` which contains ~20,982 droplets of 13 multiplexed donors,

You will most likely need to annotate your bam using ``TagReadWithGeneFunction`` (unless you have already annotated your bam with this function).
Please note that the ``\`` at the end of each line is purely for readability to put a separate parameter argument on each line.

.. admonition:: note

  If you are submitting this job to an cluster to run, you may have to bind the ``$TMPDIR`` directory used by your cluster in the singularity command.
  This will be slightly different depending on your cluster configuration.

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
.. admonition:: |:stopwatch:| Expected Resource Usage
  :class: note

  ~1.5h using a total of 5Gb memory when using 16 thread for the full :ref:`Test Dataset <TestData>` which contains ~20,982 droplets of 13 multiplexed donors,

First, we will identify the most likely singlet donor for each droplet.

.. admonition:: note

  If you are submitting this job to an cluster to run, you may have to bind the ``$TMPDIR`` directory used by your cluster in the singularity command.
  This will be slightly different depending on your cluster configuration.

Please note that the ``\`` at the end of each line is purely for readability to put a separate parameter argument on each line.

.. code-block:: bash

  singularity exec Demuxafy.sif Dropulation_AssignCellsToSamples.py --CELL_BC_FILE $BARCODES \
            --INPUT_BAM $DROPULATION_OUTDIR/possorted_genome_bam_dropulation_tag.bam \
            --OUTPUT $DROPULATION_OUTDIR/assignments.tsv.gz \
            --VCF $VCF \
            --SAMPLE_FILE $INDS \
            ${CELL_TAG:+--CELL_BARCODE_TAG $CELL_TAG} \
            ${UMI_TAG:+--MOLECULAR_BARCODE_TAG $UMI_TAG} \
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
.. admonition:: |:stopwatch:| Expected Resource Usage
  :class: note

  ~1.5h using a total of 5Gb memory when using 16 thread for the full :ref:`Test Dataset <TestData>` which contains ~20,982 droplets of 13 multiplexed donors,

Next, we will identify the likelihoods of each droplet being a doublet.

.. admonition:: note

  If you are submitting this job to an cluster to run, you may have to bind the ``$TMPDIR`` directory used by your cluster in the singularity command.
  This will be slightly different depending on your cluster configuration.

Please note that the ``\`` at the end of each line is purely for readability to put a separate parameter argument on each line.

.. code-block:: bash

  singularity exec Demuxafy.sif DetectDoublets --CELL_BC_FILE $BARCODES \
            --INPUT_BAM $DROPULATION_OUTDIR/possorted_genome_bam_dropulation_tag.bam \
            --OUTPUT $DROPULATION_OUTDIR/likelihoods.tsv.gz \
            --VCF $VCF \
            ${CELL_TAG:+--CELL_BARCODE_TAG $CELL_TAG} \
            ${UMI_TAG:+--MOLECULAR_BARCODE_TAG $UMI_TAG} \
            --SINGLE_DONOR_LIKELIHOOD_FILE $DROPULATION_OUTDIR/assignments.tsv.gz \
            --SAMPLE_FILE $INDS \
            --MAX_ERROR_RATE 0.05


If the bam annotation is successful, you will have these new files in your ``$DROPULATION_OUTDIR``:

.. code-block:: bash
  :emphasize-lines: 2,3,4

  /path/to/output/dropulation
  ├── assignments.tsv.gz
  ├── likelihoods.tsv.gz
  ├── out_vcf.vcf
  ├── out_vcf.vcf.idx
  └── possorted_genome_bam_dropulation_tag.bam


Dropulation Call
^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, we will make final assignments for each droplet based on the doublet and assignment calls.

Please note that the ``\`` at the end of each line is purely for readability to put a separate parameter argument on each line.

.. code-block:: bash

  singularity exec Demuxafy.sif dropulation_call.R --assign $DROPULATION_OUTDIR/assignments.tsv.gz \
                             --doublet $DROPULATION_OUTDIR/likelihoods.tsv.gz \
                             --out $DROPULATION_OUTDIR/updated_assignments.tsv.gz


If the bam annotation is successful, you will have these new files in your ``$DROPULATION_OUTDIR``:

.. code-block:: bash
  :emphasize-lines: 5

  /path/to/output/dropulation
  ├── assignments.tsv.gz
  ├── likelihoods.tsv.gz
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
  | 113_113         | 1327         |
  +-----------------+--------------+
  | 349_350         | 1440         |
  +-----------------+--------------+
  | 352_353         | 1562         |
  +-----------------+--------------+
  | 39_39           | 1255         |
  +-----------------+--------------+
  | 40_40           | 1082         |
  +-----------------+--------------+
  | 41_41           | 1122         |
  +-----------------+--------------+
  | 42_42           | 1365         |
  +-----------------+--------------+
  | 43_43           | 1546         |
  +-----------------+--------------+
  | 465_466         | 1084         |
  +-----------------+--------------+
  | 596_597         | 1258         |
  +-----------------+--------------+
  | 597_598         | 1515         |
  +-----------------+--------------+
  | 632_633         | 815          |
  +-----------------+--------------+
  | 633_634         | 892          |
  +-----------------+--------------+
  | 660_661         | 1364         |
  +-----------------+--------------+
  | doublet         | 3355         |
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
      | CATATGGCAGCTCGCA-1         | -44.523                | 596_597                | singlet                 | 193               | 381                     |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      | ACATACGGTCGAATCT-1         | -93.431                | 632_633                | singlet                 | 296               | 675                     |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      | GCATGCGAGATCACGG-1         | -45.708                | 41_41                  | singlet                 | 241               | 536                     |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      | CCTTACGGTAGCTCCG-1         | -21.723                | 41_41                  | singlet                 | 135               | 217                     |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      | TTTACTGCAATGAATG-1         | -26.521                | 352_353                | singlet                 | 120               | 206                     |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+
      | ...                        | ...                    | ...                    |  ...                    | ...               | ...                     |
      +----------------------------+------------------------+------------------------+-------------------------+-------------------+-------------------------+




Merging Results with Other Software Results
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.

Citation
-----------
If you used the Demuxafy platform for analysis, please reference our publication_ as well as Dropulation_.




