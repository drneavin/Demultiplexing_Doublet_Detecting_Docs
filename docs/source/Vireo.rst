.. # define a hard line break for HTML
.. |br| raw:: html

   <br />


.. _Vireo-docs:

Vireo Tutorial
===========================

.. _Vireo: https://github.com/statgen/popscle

Vireo is a flexible demultiplexing software that can demutliplex without any reference SNP genotypes, with reference SNP genotypes for a subset of the donors in the pool or no reference SNP genotypes.
If you have reference SNP genotypes for **all** of the donors in your pool, you could also use :ref:`Demuxlet <Demuxlet-docs>` or :ref:`Souporcell <Souporcell-docs>`.
If you don't have reference SNP genotypes, you could alternatively use :ref:`Freemuxlet<Freemuxlet-docs>` or :ref:`ScSplit<scSplit-docs>`.




Data
----
This is the data that you will need to have preparede to run Vireo_:


.. admonition:: Required
  :class: important

  - Common SNP genotypes vcf (``$VCF``)

    - If you have reference SNP genotypes for individuals in your pool, you can use those

      - For Vireo_ you should only have the donors that are in this pool in the vcf file

    - If you do not have reference SNP genotypes, they can be from any large population resource (i.e. 1000 Genomes or HRC)

    - Filter for common SNPs (> 5% minor allele frequency) and SNPs overlapping genes

    - Filter for common SNPs (> 5% minor allele frequency) and SNPs overlapping genes

  - Barcode file (``$BARCODES``)

  - Number of samples in pool (``$N``)
  
  - Bam file (``$BAM``)

    - Aligned single cell reads

  - Output directory (``$VIREO_OUTDIR``)
  




Run Vireo
------------
CellSNP Pileup
^^^^^^^^^^^^^^
First, you need to count the number of alleles at each SNP in each droplet using cellSNP-lite:

.. code-block:: bash

  singularity exec Demuxafy.sif cellsnp-lite -s $BAM -b $BARCODES -O $VIREO_OUTDIR -R $VCF -p 20 --minMAF 0.1 --minCOUNT 20

You can alter the ``-p``, ``--minMAF`` and ``--minCOUNT`` parameters to fit your data and your needs.
We have found these settings to work well with our data

If the pileup is successfull, you will have this new file in your ``$VIREO_OUTDIR``:

.. code-block:: bash

	.
	├── cellSNP.base.vcf
	├── cellSNP.samples.tsv
	├── cellSNP.tag.AD.mtx
	├── cellSNP.tag.DP.mtx
	└── cellSNP.tag.OTH.mtx

Additional details about outputs are available below in the :ref:`Vireo Results and Interpretation <vireo-results>`.


Vireo expects the `cellSNP.base.vcf` input to be gzipped. 
We have contacted the developers to update this issue, but until ``cellSNP-lite`` outputs, we will need to add one extra step to prepare for the demultiplexing step:

.. code-block::

	singularity exec Demuxafy.sif bgzip $VIREO_OUTDIR/cellSNP.base.vcf


Demultiplex with Vireo
^^^^^^^^^^^^^^^^^^^^^^
Next, we can use the cellSNP results to demultiplex the data with Vireo_.
As already mentioned, you can use Vireo_ with multiple different levels of reference SNP genotypes.
We've provided an example command for each of these differing amounts of donor SNP genotype data.

.. tabs::

  .. tab:: With SNP Genotype |br| Data for All Donors

    You will need to provide which genotype measure  (``$FIELD``) is provided in your donor SNP genotype file (GT, GP, or PL); default is PL.

    .. admonition:: Note

      For Vireo_ you should only have the donors that are in this pool in the vcf file.

    .. admonition:: Recommended
      :class: important

      Vireo runs more efficiently when the SNPs from the donor ``$VCF`` have been filtered for the SNPs identified by ``cellSNP-lite``.
      Therefore, it is highly recommended subset the vcf first:\.

      If your ``$VCF`` file is bgzipped (`i.e.` ends in ``.vcf.gz``), you can use ``bcftools`` to filter your ``$VCF`` for the positions in ``$VIREO_OUTDIR/cellSNP.base.vcf``:

        .. code-block::

          singularity exec Demuxafy.sif bcftools view $VCF -R $VIREO_OUTDIR/cellSNP.base.vcf.gz -Ov -o $VIREO_OUTDIR/donor_subset.vcf

      If your ``$VCF`` file is **not** bgzipped (`i.e.` ends in ``.vcf``), you can use ``bedtools`` to filter your ``$VCF`` for the positions in ``$VIREO_OUTDIR/cellSNP.base.vcf``:

        .. code-block::

          singularity exec Demuxafy.sif bedtools intersect -a $VCF -b $VIREO_OUTDIR/cellSNP.base.vcf.gz -wa -header > $VIREO_OUTDIR/donor_subset.vcf

    To run Vireo_ with reference SNP genotype data for your donors (ideally filtered as shown above):

    .. code-block::

      singularity exec Demuxafy.sif vireo -c $VIREO_OUTDIR -d $VIREO_OUTDIR/donor_subset.vcf -o $VIREO_OUTDIR -t $FIELD

  .. tab:: With SNP Genotype |br| Data for Some Donors

    .. admonition:: Note

      For Vireo_ you should only have the donors that are in this pool in the vcf file. 
      It assumes that ``$N`` is larger than the number of donors in the ``$VCF``

    .. admonition:: Recommended
      :class: important

      Vireo runs more efficiently when the SNPs from the donor ``$VCF`` have been filtered for the SNPs identified by ``cellSNP-lite``.
      Therefore, it is highly recommended subset the vcf as follows first:

        .. code-block::

          singularity exec Demuxafy.sif bcftools view $VCF -R $VIREO_OUTDIR/cellSNP.base.vcf -Oz -o $VIREO_OUTDIR/donor_subset.vcf


    .. code-block::

      singularity exec Demuxafy.sif vireo -c $VIREO_OUTDIR/cellSNPpileup.vcf.gz -d $VIREO_OUTDIR/donor_subset.vcf -o $VIREO_OUTDIR/cellSNPpileup.vcf.gz -t $FIELD -N $N

  .. tab:: Without Donor SNP |br| Genotype Data

    .. code-block::

      singularity exec Demuxafy.sif vireo -c $VIREO_OUTDIR/cellSNPpileup.vcf.gz -o $VIREO_OUTDIR -N $N

If Vireo_ is successfull, you will have these new files in your ``$VIREO_OUTDIR``:

.. code-block:: bash
  :emphasize-lines: 7,8,9,10,11,12,13

  .
  ├── cellSNP.base.vcf
  ├── cellSNP.samples.tsv
  ├── cellSNP.tag.AD.mtx
  ├── cellSNP.tag.DP.mtx
  ├── cellSNP.tag.OTH.mtx
  ├── donor_ids.tsv
  ├── donor_subset.vcf
  ├── fig_GT_distance_estimated.pdf
  ├── _log.txt
  ├── prob_doublet.tsv.gz
  ├── prob_singlet.tsv.gz
  └── summary.tsv

Additional details about outputs are available below in the :ref:`Vireo Results and Interpretation <vireo-results>`.


.. _vireo-results:

Vireo Results and Interpretation
-------------------------------------
After running the Vireo_ steps, you will have a number of files in your ``$VIREO_OUTDIR``. 
Theses are the files that most users will find the most informative:

- ``summary.tsv``

  - A summary of the droplets assigned to each donor, doublets and unassigned.

    +------------+------+
    | Var1       | Freq |
    +============+======+
    | 113_113    | 1342 |
    +------------+------+
    | 349_350    | 1475 |
    +------------+------+
    | 352_353    | 1619 |
    +------------+------+
    | 39_39      | 1309 |
    +------------+------+
    | 40_40      | 1097 |
    +------------+------+
    | 41_41      | 1144 |
    +------------+------+
    | 42_42      | 1430 |
    +------------+------+
    | 43_43      | 1561 |
    +------------+------+
    | 465_466    | 1104 |
    +------------+------+
    | 596_597    | 1271 |
    +------------+------+
    | 597_598    | 1532 |
    +------------+------+
    | 632_633    | 871  |
    +------------+------+
    | 633_634    | 967  |
    +------------+------+
    | 660_661    | 1377 |
    +------------+------+
    | doublet    | 2770 |
    +------------+------+
    | unassigned | 113  |
    +------------+------+

    - To check whether the numbe of doublets identified by Vireo_ is consistent with the expected doublet rate based on the number of droplets that you captured, you can use our `Expected Doublet Estimation Calculator <test.html>`__.


- ``donor_ids.tsv``

  - The classification of each droplet, and some droplet metrics.

    +-------------------------+---------+-----------------+-----------------+---------+--------------+------------------+
    | cell                    | donor_id|        prob_max | prob_doublet    | n_vars  | best_singlet |  best_doublet    |
    +=========================+=========+=================+=================+=========+==============+==================+
    | AAACCTGAGATAGCAT-1      | 41_41   | 1.00e+00        | 9.13e-09        | 115     | 41_41        | 40_40,41_41      |
    +-------------------------+---------+-----------------+-----------------+---------+--------------+------------------+
    | AAACCTGAGCAGCGTA-1      | 465_466 | 1.00e+00        | 5.03e-17        | 239     | 465_466      | 349_350,43_43    |
    +-------------------------+---------+-----------------+-----------------+---------+--------------+------------------+
    | AAACCTGAGCGATGAC-1      | 113_113 | 1.00e+00        | 7.57e-07        | 98      | 113_113      | 113_113,633_634  |
    +-------------------------+---------+-----------------+-----------------+---------+--------------+------------------+
    | AAACCTGAGCGTAGTG-1      | 349_350 | 1.00e+00        | 8.07e-07        | 140     | 349_350      | 349_350,597_598  |
    +-------------------------+---------+-----------------+-----------------+---------+--------------+------------------+
    | AAACCTGAGGAGTTTA-1      | 632_633 | 1.00e+00        | 5.99e-11        | 177     | 632_633      | 40_40,113_113    |
    +-------------------------+---------+-----------------+-----------------+---------+--------------+------------------+
    | AAACCTGAGGCTCATT-1      | 39_39   | 1.00e+00        | 4.44e-06        | 110     | 39_39        | 39_39,40_40      |
    +-------------------------+---------+-----------------+-----------------+---------+--------------+------------------+


Merging Results with Other Software Restults
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.

Citation
--------
If you used this workflow for analysis, please reference our paper (REFERENCE) as well as `Vireo <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1865-2>`__.