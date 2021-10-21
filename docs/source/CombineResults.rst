.. _Combine-docs:

Combining Results
=================

After you have run each of the Demultiplexing and Doublet Detecting softwares you would like, it is helpful to convert them to similar nomenclarture and combine the results into a single dataframe.
In addition, we have found it helpful to generate summaries of each of the combinations of softwares identified.
To help streamline this process, we have provided a script that will easily integrate all the softwares you have run into a single dataframe and to generate a summary file for all the software combinations and if you ran demultiplexing softwares, it will also generate a demultiplexing summary file for the individual and cluster assignments from the demultiplexing softwares.


Data
-----
In order to use our script to combine the results from the various demultiplexing and doublet detecting softwares, you need the following:

.. admonition:: Required
  :class: important

  - Output directory (``$OUTDIR``)

  - Path to results of each of the softwares you would like to merge into a single dataframe.

    - You need to provide the path to at least one software result, otherwise, it will not run.


Merging Results with Combine_Results.r
--------------------------------------

The script has multiple options to provide the paths to each of the software results you would like to run.
To see each of the options, simply run:

.. code-block:: bash

  singularity exec Demuxafy.sif Combine_Results.R -h

Providing the possible parameter options:

.. code-block:: bash

  usage: Combine_Results.R
       [-h] -o OUT [-d DEMUXLET] [-f FREEMUXLET] [-g FREEMUXLET_ASSIGNMENTS]
       [-s SCSPLIT] [-w SCSPLIT_ASSIGNMENTS] [-u SOUPORCELL]
       [-x SOUPORCELL_ASSIGNMENTS] [-v VIREO] [-e DOUBLETDECON]
       [-t DOUBLETDETECTION] [-i DOUBLETFINDER] [-n SCDBLFINDER] [-c SCDS]
       [-r SCRUBLET] [-l SOLO]

  optional arguments:
    -h, --help            show this help message and exit
    -o OUT, --out OUT     The file where results will be saved
    -d DEMUXLET, --demuxlet DEMUXLET
                          Path to demuxlet results. Only use this option if you
                          want to include the demuxlet results.
    -f FREEMUXLET, --freemuxlet FREEMUXLET
                          Path to freemuxlet results. Only use this option if
                          you want to include the freemuxlet results.
    -g FREEMUXLET_ASSIGNMENTS, --freemuxlet_assignments FREEMUXLET_ASSIGNMENTS
                          Path to freemuxlet cluster-to-individual assignments.
                          Only use this option if have used reference SNP
                          genotypes to assign individuals to clusters for the
                          freemuxlet results.
    -s SCSPLIT, --scSplit SCSPLIT
                          Path to scSplit results. Only use this option if you
                          want to include the scSplit results.
    -w SCSPLIT_ASSIGNMENTS, --scSplit_assignments SCSPLIT_ASSIGNMENTS
                          Path to scSplit cluster-to-individual assignments.
                          Only use this option if you have used reference SNP
                          genotypes to assign individuals to clusters for the
                          scSplit results.
    -u SOUPORCELL, --souporcell SOUPORCELL
                          Path to souporcell results. Only use this option if
                          you want to include the souporcell results.
    -x SOUPORCELL_ASSIGNMENTS, --souporcell_assignments SOUPORCELL_ASSIGNMENTS
                          Path to souporcell cluster-to-individual assignments.
                          Only use this option if you have used reference SNP
                          genotypes to assign individuals to clusters for the
                          souporcell results.
    -v VIREO, --vireo VIREO
                          Path to vireo results. Only use this option if you
                          want to include the vireo results.
    -e DOUBLETDECON, --DoubletDecon DOUBLETDECON
                          Path to DoubletDecon results. Only use this option if
                          you want to include the DoubletDecon results.
    -t DOUBLETDETECTION, --DoubletDetection DOUBLETDETECTION
                          Path to DoubletDetection results. Only use this option
                          if you want to include the DoubletDetection results.
    -i DOUBLETFINDER, --DoubletFinder DOUBLETFINDER
                          Path to DoubletFinder results. Only use this option if
                          you want to include the DoubletFinder results.
    -n SCDBLFINDER, --scDblFinder SCDBLFINDER
                          Path to scDblFinder results. Only use this option if
                          you want to include the scDblFinder results.
    -c SCDS, --scds SCDS  Path to scds results. Only use this option if you want
                          to include the scds results.
    -r SCRUBLET, --scrublet SCRUBLET
                          Path to scrublet results. Only use this option if you
                          want to include the scrublet results.
    -l SOLO, --solo SOLO  Path to solo results. Only use this option if you want
                          to include the solo results.

An example command that combines :ref:`Demuxlet <Demuxlet-docs>` results, :ref:`Souporcell <Souporcell-docs>` results, :ref:`Solo <Solo-docs>` results and :ref:`Scds <scds-docs>` results would look like this:

.. code-block:: bash

  singularity exec Demuxafy.sif Combine_Results.R \
    -o $OUTDIR/combined_results.tsv \
    --demuxlet $DEMUXLET_OUTDIR \
    --souporcell $DEMUXLET_OUTDIR \
    --solo $SOLO_OUTDIR \
    --scds $SCDS_OUTDIR \


.. admonition:: Note

  The path to the directories will work if the file names are the expected file names.
  However, if you used a different file naming convention or changed the names, you can also provide the full path to the exact file for each software.


Results and Interpretation
--------------------------
After running the ``Combine_Results.R`` script, you should have three (or two if you didn't have any demultiplexing softwares)

.. code-block:: bash

  .
  ├── combined_results_demultiplexing_summary.tsv
  ├── combined_results_summary.tsv
  └── combined_results.tsv

.. admonition:: Note

  You will only have the ``combined_results_demultiplexing_summary.tsv`` file if you included demultiplexing softwares.

Here's a deeper look at the contents of each of these cells:

 - ``combined_results.tsv``
  
    - Has the selected results combined; only including key columns.

      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+
      | Barcode            |Demuxlet_DropletType | Demuxlet_Individual_Assignment | Souporcell_Individual_Assignment | Souporcell_Cluster | Souporcell_DropletType | scds_score         | scds_DropletType | solo_DropletType | solo_DropletScore |
      +====================+=====================+================================+==================================+====================+========================+====================+==================+==================+===================+
      | AAACCTGAGATAGCAT-1 | singlet             | 41_41                          | 41_41                            | 6                  | singlet                | 0.116344358493288  | singlet          | singlet          | -8.442187         |
      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+
      | AAACCTGAGCAGCGTA-1 | singlet             | 465_466                        | 465_466                          | 11                 | singlet                | 0.539856378453988  | singlet          | singlet          | -2.8096201        |
      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+
      | AAACCTGAGCGATGAC-1 | singlet             | 113_113                        | 113_113                          | 5                  | singlet                | 0.0237184380134577 | singlet          | singlet          | -2.8949203        |
      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+
      | AAACCTGAGCGTAGTG-1 | singlet             | 349_350                        | 349_350                          | 3                  | singlet                | 0.163695865366576  | singlet          | singlet          | -5.928284         |
      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+
      | AAACCTGAGGAGTTTA-1 | singlet             | 632_633                        | 632_633                          | 7                  | singlet                | 0.11591462421927   | singlet          | doublet          | 0.2749935         |
      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+
      | AAACCTGAGGCTCATT-1 | singlet             | 39_39                          | 39_39                            | 12                 | singlet                | 0.0479944175570073 | singlet          | singlet          | -5.2726507        |
      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+
      | AAACCTGAGGGCACTA-1 | singlet             | 465_466                        | 465_466                          | 11                 | singlet                | 0.374426050641161  | singlet          | singlet          | -0.65760195       |
      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+
      | AAACCTGAGTAATCCC-1 | singlet             | 660_661                        | 660_661                          | 4                  | singlet                | 0.247842972104563  | singlet          | singlet          | -3.5948637        |
      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+
      | AAACCTGAGTAGCCGA-1 | doublet             | doublet                        | unassigned                       | unassigned         | unassigned             | 0.342998285281922  | singlet          | singlet          | -0.50507957       |
      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+
      | ...                | ...                 | ...                            | ...                              | ...                | ...                    | ...                | ...              | ...              | ...               |
      +--------------------+---------------------+--------------------------------+----------------------------------+--------------------+------------------------+--------------------+------------------+------------------+-------------------+

  - ``combined_results_summary.tsv``

    - The number of each of the combinations of the software cell type classifications

    +----------------------+-------------------------+-------------------+-------------------+-------+
    | Demuxlet_DropletType | Souporcell_DropletType  | scds_DropletType  | solo_DropletType  | N     |
    +======================+=========================+===================+===================+=======+
    | singlet              | singlet                 | singlet           | singlet           | 16193 |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | doublet                 | doublet           | doublet           | 1714  |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | singlet              | singlet                 | singlet           | doublet           | 947   |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | doublet                 | singlet           | singlet           | 468   |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | singlet              | singlet                 | doublet           | singlet           | 392   |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | singlet              | singlet                 | doublet           | doublet           | 345   |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | doublet                 | singlet           | doublet           | 335   |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | singlet                 | singlet           | singlet           | 171   |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | doublet                 | doublet           | singlet           | 169   |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | singlet                 | doublet           | doublet           | 114   |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | singlet                 | singlet           | doublet           | 44    |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | singlet                 | doublet           | singlet           | 18    |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | singlet              | doublet                 | singlet           | singlet           | 17    |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | singlet              | unassigned              | singlet           | singlet           | 13    |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | unassigned              | singlet           | singlet           | 11    |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | singlet              | doublet                 | doublet           | doublet           | 9     |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | singlet              | doublet                 | singlet           | doublet           | 6     |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | singlet              | doublet                 | doublet           | singlet           | 5     |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | unassigned              | singlet           | doublet           | 4     |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | unassigned              | doublet           | doublet           | 3     |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | doublet              | unassigned              | doublet           | singlet           | 2     |
    +----------------------+-------------------------+-------------------+-------------------+-------+
    | unassigned           | unassigned              | singlet           | singlet           | 2     |
    +----------------------+-------------------------+-------------------+-------------------+-------+

  - ``combined_results_demultiplexing_summary.tsv``

    - Summary of the number of each of the combination of classifications by demultiplexing software:

      +---------------------------------+-----------------------------------------+------+
      | Demuxlet_Individual_Assignment  | Souporcell_Individual_Assignment        | N    |
      +=================================+=========================================+======+
      | doublet                         | doublet                                 | 2686 |
      +---------------------------------+-----------------------------------------+------+
      | 352_353                         | 352_353                                 | 1603 |
      +---------------------------------+-----------------------------------------+------+
      | 43_43                           | 43_43                                   | 1547 |
      +---------------------------------+-----------------------------------------+------+
      | 597_598                         | 597_598                                 | 1510 |
      +---------------------------------+-----------------------------------------+------+
      | 349_350                         | 349_350                                 | 1450 |
      +---------------------------------+-----------------------------------------+------+
      | 42_42                           | 42_42                                   | 1417 |
      +---------------------------------+-----------------------------------------+------+
      | 660_661                         | 660_661                                 | 1358 |
      +---------------------------------+-----------------------------------------+------+
      | 113_113                         | 113_113                                 | 1333 |
      +---------------------------------+-----------------------------------------+------+
      | 39_39                           | 39_39                                   | 1289 |
      +---------------------------------+-----------------------------------------+------+
      | ...                             | ...                                     | ...  |
      +---------------------------------+-----------------------------------------+------+



Citation
--------
If you used this workflow for analysis, please reference our paper (REFERENCE).