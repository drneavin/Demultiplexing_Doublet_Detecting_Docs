.. _Demuxalot-docs:


Demuxalot
===========================

.. _Demuxalot: https://pypi.org/project/demuxalot/
.. _preprint: https://www.biorxiv.org/content/10.1101/2022.03.07.483367v1




Demuxalot_ is a genotype demultiplexing software that requires reference genotypes to be available for each individual in the pool. 
Therefore, if you don't have reference genotypes, you may want to demultiplex with one of the softwares that do not require reference genotype data
(:ref:`Freemuxlet <Freemuxlet-docs>`, :ref:`scSplit <scSplit-docs>`, :ref:`Souporcell <Souporcell-docs>` or :ref:`Vireo <Vireo-docs>`)


Data
----
This is the data that you will need to have prepare to run Demuxalot_:

.. admonition:: Required
  :class: important

  - Reference SNP genotypes for each individual (``$VCF``)

    - Filter for common SNPs (> 5% minor allele frequency) and SNPs overlapping genes

  - Barcode file (``$BARCODES``)

  - Bam file (``$BAM``)

    - Aligned single cell reads

  - Output directory (``$DEMUXALOT_OUTDIR``)

  - A text file with the individual ids (``$INDS``)
    
      - File containing the individual ids (separated by line) as they appear in the vcf file

      - For example, this is the :download:`individual file <_download_files/Individuals.txt>` for our example dataset



Run Demuxalot
----------------

.. admonition:: Example Variable Settings
  :class: grey

    Below is an example of the variables that we can set up to be used in the command below.
    These are files provided as a :ref:`test dataset <TestData>` available in the :ref:`Data Preparation Documentation <DataPrep-docs>`
    Please replace ``/path/to`` with the full path to your data directory.

    .. code-block:: bash

        VCF=/path/to/TestData4PipelineFull/test_dataset.vcf
        BARCODES=/path/to/TestData4PipelineFull/test_dataset/outs/filtered_gene_bc_matrices/Homo_sapiens_GRCh38p10/barcodes.tsv
        BAM=/path/to/test_dataset/possorted_genome_bam.bam
        DEMUXALOT_OUTDIR=/path/to/output/demuxalot
        INDS=/path/to/TestData4PipelineFull/donor_list.txt

Demuxalot_ can be run with refinement if desired - this means that there is another step that is run to refine the genotypes from the data which makes the results slightly more accurate.
You can choose to run Demuxalot_ with or without refinement:

Demultiplex with Demuxalot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. admonition:: |:stopwatch:| Expected Resource Usage
  :class: note

  ~2.5h using a total of 81Gb memory when using 32 threads for the full :ref:`Test Dataset <TestData>` which contains ~20,982 droplets of 13 multiplexed donors,


.. tabs::

  .. tab:: With Refinement

    This will run the first phase of Demuxalot_ as well as the subsequent refinement:

    .. code-block:: bash

      singularity exec Demuxafy.sif Demuxalot.py \
              -b $BARCODES \
              -a $BAM \
              -n $INDS \
              -v $VCF \
              -o $DEMUXALOT_OUTDIR \
              -r True

    .. admonition:: HELP! It says my file/directory doesn't exist!
      :class: dropdown

      If you receive an error indicating that a file or directory doesn't exist but you are sure that it does, this is likely an issue arising from Singularity.
      This is easy to fix.
      The issue and solution are explained in detail in the :ref:`Notes About Singularity Images <Singularity-docs>`

    If Demuxalot_ is successful, you will have these new files in your ``$DEMUXALOT_OUTDIR``:

    .. code-block:: bash

      /path/to/output/demuxalot
      ├── assignments_refined.tsv.gz
      ├── assignments.tsv.gz
      ├── likelihoods_refined.tsv.gz
      ├── likelihoods.tsv.gz
      ├── posterior_probabilities_refined.tsv.gz
      └── posterior_probabilities.tsv.gz


  .. tab:: Without Refinement

    This will run the first phase of Demuxalot_ only without any refinement:

    .. code-block:: bash

      singularity exec Demuxafy.sif Demuxalot.py \
              -b $BARCODES \
              -a $BAM \
              -n $INDS \
              -v $VCF \
              -o $DEMUXALOT_OUTDIR \
              -r False

    .. admonition:: HELP! It says my file/directory doesn't exist!
      :class: dropdown

      If you receive an error indicating that a file or directory doesn't exist but you are sure that it does, this is likely an issue arising from Singularity.
      This is easy to fix.
      The issue and solution are explained in detail in the :ref:`Notes About Singularity Images <Singularity-docs>`

    If Demuxalot_ is successful, you will have these new files in your ``$DEMUXALOT_OUTDIR``:

    .. code-block:: bash

      /path/to/output/demuxalot
      ├── assignments.tsv.gz
      ├── likelihoods.tsv.gz
      └── posterior_probabilities.tsv.gz

Additional details about outputs are available below in the :ref:`Demuxalot Results and Interpretation <demuxalot-results>`.


Demuxalot Summary
^^^^^^^^^^^^^^^^^^^
We have provided a script that will summarize the number of droplets classified as doublets, ambiguous and assigned to each donor by Demuxalot_ and write it to the ``$DEMUXALOT_OUTDIR``. 
You can run this to get a fast and easy summary of your results by providing the path to your result file:

.. tabs::

  .. tab:: With Refinement

    .. code-block:: bash

      singularity exec Demuxafy.sif bash demuxalot_summary.sh $DEMUXALOT_OUTDIR/assignments_refined.tsv.gz


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

      singularity exec Demuxafy.sif bash Demuxalot_summary.sh $DEMUXLET_OUTDIR/assignments_refined.tsv.gz > $DEMUXLET_OUTDIR/demuxalot_summary.tsv


  .. tab:: Without Refinement

    .. code-block:: bash

      singularity exec Demuxafy.sif bash Demuxalot_summary.sh $DEMUXALOT_OUTDIR/assignments.tsv.gz


    which will return:

      +-----------------+--------------+
      | Classification  | Assignment N |
      +=================+==============+
      | 113_113         | 1344         |
      +-----------------+--------------+
      | 349_350         | 1463         |
      +-----------------+--------------+
      | 352_353         | 1619         |
      +-----------------+--------------+
      | 39_39           | 1306         |
      +-----------------+--------------+
      | 40_40           | 1082         |
      +-----------------+--------------+
      | 41_41           | 1129         |
      +-----------------+--------------+
      | 42_42           | 1437         |
      +-----------------+--------------+
      | 43_43           | 1553         |
      +-----------------+--------------+
      | 465_466         | 1091         |
      +-----------------+--------------+
      | 596_597         | 1267         |
      +-----------------+--------------+
      | 597_598         | 1523         |
      +-----------------+--------------+
      | 632_633         | 872          |
      +-----------------+--------------+
      | 633_634         | 961          |
      +-----------------+--------------+
      | 660_661         | 1371         |
      +-----------------+--------------+
      | doublet         | 2964         |
      +-----------------+--------------+




    or you can write it straight to a file:

    .. code-block:: bash

      singularity exec Demuxafy.sif bash Demuxalot_summary.sh $DEMUXLET_OUTDIR/assignments.tsv.gz > $DEMUXLET_OUTDIR/demuxalot_summary.tsv



.. admonition:: Note

  To check if these numbers are consistent with the expected doublet rate in your dataset, you can use our `Doublet Estimation Calculator <test.html>`__.


.. _demuxalot-results:

Demuxalot Results and Interpretation
-----------------------------------------
After running the Demuxalot_ steps and summarizing the results, you will have a number of files from some of the intermediary steps. 
These are the files that most users will find the most informative:

  - ``assignments.tsv.gz`` (and ``assignments_refined.tsv.gz`` if you indicated ``-r True``)

    - The droplet assignment for each barcode:

      +---------------------------+------------------+
      | BARCODE                   | 0                |
      +===========================+==================+
      | AAACCTGAGATAGCAT-1        | 41_41            |
      +---------------------------+------------------+
      | AAACCTGAGCAGCGTA-1        | 465_466          |
      +---------------------------+------------------+
      | AAACCTGAGCGATGAC-1        | 113_113          |
      +---------------------------+------------------+
      | AAACCTGAGCGTAGTG-1        | 349_350          |
      +---------------------------+------------------+
      | AAACCTGAGGAGTTTA-1        | 632_633          |
      +---------------------------+------------------+
      | AAACCTGAGGCTCATT-1        | 39_39            |
      +---------------------------+------------------+
      | ...                       | ...              |
      +---------------------------+------------------+


  - ``likelihoods_refined.tsv.gz`` or ``likelihoods.tsv.gz``:

    - The likelihood probabilities for each donor and doublet combination for each droplet

  - ``posterior_probabilities_refined.tsv.gz`` or ``posterior_probabilities_refined.tsv.gz``

    - The posterior probabilities for each donor or doublet combination for each droplet


Merging Results with Other Software Results
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.

Citation
--------
If you used the Demuxafy platform for analysis, please reference our preprint_ as well as Demuxalot_.


