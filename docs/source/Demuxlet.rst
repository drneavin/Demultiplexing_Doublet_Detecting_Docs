.. _Demuxlet-docs:

Demuxlet Tutorial
===========================

.. _Demuxlet: https://github.com/statgen/popscle

Demuxlet_ is a genotype demultiplexing software that requires reference genotypes to be available for each individual in the pool. 
Therefore, if you don't have reference genotypes, you may want to demultiplex with one of the softwares that do not require reference genotype data
(:ref:`Freemuxlet <Freemuxlet-docs>`, :ref:`scSplit <scSplit-docs>`, :ref:`Souporcell <Souporcell-docs>` or :ref:`Vireo<Vireo-docs>`)


Data
----
This is the data that you will need to have preparede to run Demuxlet_:

.. admonition:: Required
  :class: important

  - Reference SNP genotypes for each individual (``$VCF``)

    - Filter for common SNPs (> 5% minor allele frequency) and SNPs overlapping genes

  - Genotype field in ``$VCF`` (``$FIELD``)

    - This is ``GP`` by default but could also be ``GT`` others

  - Barcode file (``$BARCODES``)

  - Bam file (``$BAM``)

    - Aligned single cell reads

  - Output directory (``$DEMUXLET_OUTDIR``)


.. admonition:: Optional

    - A text file with the individual ids (``$INDS``)
    
      - File containing the individual ids (separated by line) as they appear in the vcf file

      - For example, this is the :download:`individual file <_download_files/Individuals.txt>` for our example dataset


Run Demuxlet
------------
Poscle Pileup
^^^^^^^^^^^^^
First we will need to identify the number of reads from each allele at each SNP location:

.. code-block:: bash

  singularity exec Demuxafy.sif popscle dsc-pileup --sam $BAM --vcf $VCF --group-list $BARCODES --out $DEMUXLET_OUTDIR/popscle

If the pileup is successfull, you will have these files in your ``$DEMUXLET_OUTDIR``:

.. code-block:: bash

  .
  ├── pileup.cel.gz
  ├── pileup.plp.gz
  ├── pileup.umi.gz
  └── pileup.var.gz

Additional details about outputs are available below in the :ref:`Demuxlet Results and Interpretation <demuxlet-results>`.


Popscle Demuxlet
^^^^^^^^^^^^^^^^
Once you have run ``popscle pileup``, you can demultiplex your samples:

.. tabs::

  .. tab:: With ``$INDS`` file
    
    The ``$INDS`` file allows demuxlet to only consider the individual in this pool

    .. code-block:: bash

      singularity exec Demuxafy.sif popscle demuxlet --plp $DEMUXLET_OUTDIR/pileup --vcf $VCF --field $FIELD --group-list $BARCODES --geno-error-coeff 1.0 --geno-error-offset 0.05 --out $DEMUXLET_OUTDIR/demuxlet --sm-list $INDS

  .. tab:: Without ``$INDS`` file

    This will use all the individuals in your reference SNP genotype ``$VCF``. 
    If your ``$VCF`` only has the individuals multiplexed in your pool, then the ``$INDS`` file is not required.

    .. code-block:: bash

      singularity exec Demuxafy.sif popscle demuxlet --plp $DEMUXLET_OUTDIR/pileup --vcf $VCF --field $FIELD --group-list $BARCODES --geno-error-coeff 1.0 --geno-error-offset 0.05 --out $DEMUXLET_OUTDIR/demuxlet

If demuxlet is successfull, you will have these new files in your ``$DEMUXLET_OUTDIR``:

.. code-block:: bash
  :emphasize-lines: 2

  .
  ├── demuxlet.best
  ├── pileup.cel.gz
  ├── pileup.plp.gz
  ├── pileup.umi.gz
  └── pileup.var.gz

Additional details about outputs are available below in the :ref:`Demuxlet Results and Interpretation <demuxlet-results>`.


Demuxlet Summary
^^^^^^^^^^^^^^^^
We have provided a script that will summarize the number of droplets classified as doublets, ambiguous and assigned to each donor by Demuxlet_ and write it to the ``$DEMUXLET_OUTDIR``. 
You can run this to get a fast and easy summary of your results with:

.. code-block:: bash

  singularity exec Demuxafy.sif bash Demuxlet_summary.sh $DEMUXLET_OUTDIR

If the demuxlet summary is successfull, you will have this new file in your ``$DEMUXLET_OUTDIR``:

.. code-block:: bash
  :emphasize-lines: 3

  .
  ├── demuxlet.best
  ├── Demuxlet_summary.tsv
  ├── pileup.cel.gz
  ├── pileup.plp.gz
  ├── pileup.umi.gz
  └── pileup.var.gz

Additional details about outputs are available below in the :ref:`Demuxlet Results and Interpretation <demuxlet-results>`.




.. _demuxlet-results:

Demuxlet Results and Interpretation
-----------------------------------
After running the Demuxlet_ steps and summarizing the results, you will have a number of files from some of the intermediary steps. 
Theses are the files that most users will find the most informative:

  - ``Demuxlet_summary.tsv``

    - Summary of the droplets asignmened to each donor, doublets or unassigned.
     
      +-----------------+--------------+
      | Classification  | Assignment N |
      +=================+==============+
      | 349_350         | 1458         |
      +-----------------+--------------+
      | 113_113         | 1334         |
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
      | 633_634         | 960          |
      +-----------------+--------------+
      | 632_633         | 868          |
      +-----------------+--------------+
      | 660_661         | 1362         |
      +-----------------+--------------+
      | doublet         | 3053         |
      +-----------------+--------------+

    - To check whether the number of predicted doublets is in line with the expected number, you can use our `Doublet Estimation Calculator <test.html>`__

  - ``demuxlet.best``

    - Metrics for each droplet including the singelt, doublet or ambiguous assignment (``DROPLET.TYPE``), final assignment (``BEST.GUESS``), log likelihood of the final assignment (``BEST.LLK``) and other QC metrics.

      +---------+--------------------+----------+-----------+--------------+-------------------------+---------+-------------------------+---------+--------------------+----------------+---------------+---------------+--------------+---------------+---------------+-------------------------+-------------------------+----------------+-------------------+
      | INT_ID  | BARCODE            | NUM.SNPS | NUM.READS | DROPLET.TYPE | BEST.GUESS              |BEST.LLK |       NEXT.GUESS        |NEXT.LLK | DIFF.LLK.BEST.NEXT | BEST.POSTERIOR | SNG.POSTERIOR | SNG.BEST.GUESS| SNG.BEST.LLK | SNG.NEXT.GUESS| SNG.NEXT.LLK  | SNG.ONLY.POSTERIOR      | DBL.BEST.GUESS          |  DBL.BEST.LLK  |  DIFF.LLK.SNG.DBL |
      +=========+====================+==========+===========+==============+=========================+=========+=========================+=========+====================+================+===============+===============+==============+===============+===============+=========================+=========================+================+===================+
      | 0       | AAACCTGAGATAGCAT-1 |      170 |     231   |     SNG      | 41_41,41_41,0.00        | -29.42  | 40_40,41_41,0.50        | -39.12  | 9.70               | -33            |   1           | 41_41         | -29.42       |  597_598      | -76.24        | 0.00000                 | 40_40,41_41,0.50        | -39.12         | 9.70              |
      +---------+--------------------+----------+-----------+--------------+-------------------------+---------+-------------------------+---------+--------------------+----------------+---------------+---------------+--------------+---------------+---------------+-------------------------+-------------------------+----------------+-------------------+
      | 1       | AAACCTGAGCAGCGTA-1 |      325 |     583   |     SNG      | 465_466,465_466,0.00    | -70.61  | 42_42,465_466,0.50      | -94.85  | 24.24              | -74            |   1           | 465_466       | -70.61       |  42_42        | -166.61       | 0.00000                 | 42_42,465_466,0.50      | -94.85         | 24.24             |
      +---------+--------------------+----------+-----------+--------------+-------------------------+---------+-------------------------+---------+--------------------+----------------+---------------+---------------+--------------+---------------+---------------+-------------------------+-------------------------+----------------+-------------------+
      | 2       | AAACCTGAGCGATGAC-1 |      147 |     227   |     SNG      | 113_113,113_113,0.00    | -25.05  | 39_39,113_113,0.50      | -29.85  | 4.80               | -28            |   1           | 113_113       | -25.05       |  349_350      | -51.63        | 0.00000                 | 39_39,113_113,0.50      | -29.85         | 4.80              |
      +---------+--------------------+----------+-----------+--------------+-------------------------+---------+-------------------------+---------+--------------------+----------------+---------------+---------------+--------------+---------------+---------------+-------------------------+-------------------------+----------------+-------------------+
      | 3       | AAACCTGAGCGTAGTG-1 |      180 |     235   |     SNG      | 349_350,349_350,0.00    | -33.14  | 349_350,632_633,0.50    | -44.78  | 11.64              | -36            |   1           | 349_350       | -33.14       |  632_633      | -77.41        | 0.00000                 | 349_350,632_633,0.50    | -44.78         | 11.64             |
      +---------+--------------------+----------+-----------+--------------+-------------------------+---------+-------------------------+---------+--------------------+----------------+---------------+---------------+--------------+---------------+---------------+-------------------------+-------------------------+----------------+-------------------+
      | 4       | AAACCTGAGGAGTTTA-1 |      248 |     444   |     SNG      | 632_633,632_633,0.00    | -54.79  | 352_353,632_633,0.50    | -72.23  | 17.43              | -58            |   1           | 632_633       | -54.79       |  633_634      | -163.24       | 0.00000                 | 352_353,632_633,0.50    | -72.23         | 17.43             |
      +---------+--------------------+----------+-----------+--------------+-------------------------+---------+-------------------------+---------+--------------------+----------------+---------------+---------------+--------------+---------------+---------------+-------------------------+-------------------------+----------------+-------------------+
      | ...     | ...                | ...      | ...       | ...          | ...                     | ...     | ...                     | ...     | ...                |  ...           | ...           | ...           | ...          | ...           | ...           | ...                     | ...                     | ...            | ...               |
      +---------+--------------------+----------+-----------+--------------+-------------------------+---------+-------------------------+---------+--------------------+----------------+---------------+---------------+--------------+---------------+---------------+-------------------------+-------------------------+----------------+-------------------+


Merging Results with Other Software Restults
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.

Citation
--------
If you used this workflow for analysis, please reference our paper (REFERENCE) as well as `Demuxlet <https://www.nature.com/articles/nbt.4042>`__.