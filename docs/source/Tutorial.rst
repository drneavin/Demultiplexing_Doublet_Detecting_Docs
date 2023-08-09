Tutorial
==========

In this tutorial, we'll analyze the :ref:`Test Dataset <TestData>` which contains 13 multiplexed donors and ~20,982 droplets.

This tutorial will take you through the typical steps for demultiplexing data with Demuxafy. 
The process would be the same for non-multiplexed data for doublet detection but would be a different combination of softwares.
1. Select appropriate software combination
1. Run each selected software
1. Combine the results and call a final assignment for each droplet


Selecting Software Combination
-------------------------------
First, we'll identify the softwares we should run using the :ref:`Software Selection Tool <SoftwareSelection-docs>`.

Since we have a genetically multiplexed sample with reference SNP genotypes for each sample in the pool and less than 16 donors, we will use the following recommended softwares:

  - :ref:`Demuxlet <Demuxlet-docs>`
  - :ref:`Souporcell <Souporcell-docs>`
  - :ref:`Vireo <Vireo-docs>`
  - :ref:`Scds <scds-docs>`




Run Each Software
-------------------------------

Demuxlet
^^^^^^^^^


Souporcell
^^^^^^^^^



Vireo
^^^^^^^^^



Scds
^^^^^^^^^





Merging Results with Combine_Results.R
--------------------------------------

An example command that combines :ref:`Demuxlet <Demuxlet-docs>` results, :ref:`Souporcell <Souporcell-docs>` results, :ref:`Solo <Solo-docs>` results and :ref:`Scds <scds-docs>` results would look like this:
There are a two different options for using this script:
1. Combining the results **and** calling the droplet type through the combination of the softwares.
1. Combining the results without any joint droplet calling. You  might choose this if you just want to see how the different softwares perform on your data before deciding which to move forward with for final joint droplet calling.

First, let's assign the variables that will be used to execute each step.

.. admonition:: Example Variable Settings
  :class: grey

    Below is an example of the variables that we can set up to be used in the command below.
    These are files provided as a :ref:`test dataset <TestData>` available in the :ref:`Data Preparation Documentation <DataPrep-docs>`
    Please replace paths with the full path to data on your system.

    .. code-block:: bash

      OUTDIR=/path/to/output/combined
      DEMUXLET_OUTDIR=/path/to/output/demuxlet
      SOUPORCELL_OUTDIR=/path/to/output/souporcell
      VIREO_OUTDIR=/path/to/output/vireo
      SCDS_OUTDIR=/path/to/output/scds


.. tabs::

  .. tab:: Combine Results + Joint Droplet Calls

    The first option is to select a method to make joint calls on the individual assignment and singlet-doublet droplet types using the softwares included.

    .. code-block:: bash

      singularity exec Demuxafy.sif Combine_Results.R \
        -o $OUTDIR/combined_results.tsv \
        --demuxlet $DEMUXLET_OUTDIR \
        --souporcell $SOUPORCELL_OUTDIR \
        --solo $SOLO_OUTDIR \
        --scds $SCDS_OUTDIR \
        --method "MajoritySinglet" ## there are other methods that can also be used, please see the help message above for the other options





Results and Interpretation
--------------------------
After running the ``Combine_Results.R`` script, you should have two, three or four files depending on if you used demultiplexing softwares and if you used joint droplet calling.
Here, we show the results for the above example that also provides combined calls with the "MajoritySinglet" calls.

.. code-block:: bash

  /path/to/output/combined
  ├── combined_results_demultiplexing_summary.tsv
  ├── combined_resultsSinglets_upset.pdf
  ├── combined_results_summary.tsv
  ├── combined_results.tsv
  └── combined_results_w_combined_assignments.tsv
  

Here's a deeper look at the contents of each of these results:

  - ``combined_resultsSinglets_upset.pdf``

    - This is an upset figure of the droplets which are colored by their finall individual or doublet classification.

    - A filled circle indicates the that those droplets are classified as singlets by that method while empty circles indicate a doublet classification by that software

    .. image:: _figures/combined_resultsSinglets_upset.png

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

      +--------------------------------+-----------------------------------------+------+
      |Demuxlet_Individual_Assignment  | Souporcell_Individual_Assignment        | N    |
      +================================+=========================================+======+
      |doublet                         | doublet                                 | 2706 |
      +--------------------------------+-----------------------------------------+------+
      |352_353                         | 352_353                                 | 1603 |
      +--------------------------------+-----------------------------------------+------+
      |43_43                           | 43_43                                   | 1547 |
      +--------------------------------+-----------------------------------------+------+
      |597_598                         | 597_598                                 | 1510 |
      +--------------------------------+-----------------------------------------+------+
      |349_350                         | 349_350                                 | 1450 |
      +--------------------------------+-----------------------------------------+------+
      |42_42                           | 42_42                                   | 1417 |
      +--------------------------------+-----------------------------------------+------+
      |660_661                         | 660_661                                 | 1358 |
      +--------------------------------+-----------------------------------------+------+
      |113_113                         | 113_113                                 | 1333 |
      +--------------------------------+-----------------------------------------+------+
      |39_39                           | 39_39                                   | 1289 |
      +--------------------------------+-----------------------------------------+------+
      |...                             | ...                                     | ...  |
      +--------------------------------+-----------------------------------------+------+

  - ``combined_results_w_combined_assignments.tsv``

    - Dataframe combining all the software results together + combined assignment based on selected method:

    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+
    | Barcode                 | Demuxlet_DropletType    | Demuxlet_Individual_Assignment  | Souporcell_Cluster      | Souporcell_Individual_Assignment  | Souporcell_DropletType  | scds_score            | scds_DropletType        | solo_DropletType        | solo_DropletScore | MajoritySinglet_DropletType     | MajoritySinglet_Individual_Assignment|
    +=========================+=========================+=================================+=========================+===================================+=========================+=======================+=========================+=========================+===================+=================================+======================================+
    | AAACCTGAGATAGCAT-1      | singlet                 | 41_41                           | 6                       | 41_41                             | singlet                 | 0.116344358493288     | singlet                 | singlet                 | -8.442187         | singlet                         |  41_41                               |
    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+
    | AAACCTGAGCAGCGTA-1      | singlet                 | 465_466                         | 11                      | 465_466                           | singlet                 | 0.539856378453988     | singlet                 | singlet                 | -2.8096201        | singlet                         |  465_466                             |
    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+
    | AAACCTGAGCGATGAC-1      | singlet                 | 113_113                         | 5                       | 113_113                           | singlet                 | 0.0237184380134577    | singlet                 | singlet                 | -2.8949203        | singlet                         |  113_113                             |
    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+
    | AAACCTGAGCGTAGTG-1      | singlet                 | 349_350                         | 3                       | 349_350                           | singlet                 | 0.163695865366576     | singlet                 | singlet                 | -5.928284         | singlet                         |  349_350                             |
    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+
    | AAACCTGAGGAGTTTA-1      | singlet                 | 632_633                         | 7                       | 632_633                           | singlet                 | 0.11591462421927      | singlet                 | doublet                 | 0.2749935         | singlet                         |  632_633                             |
    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+
    | AAACCTGAGGCTCATT-1      | singlet                 | 39_39                           | 12                      | 39_39                             | singlet                 | 0.0479944175570073    | singlet                 | singlet                 | -5.2726507        | singlet                         |  39_39                               |
    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+
    | AAACCTGAGGGCACTA-1      | singlet                 | 465_466                         | 11                      | 465_466                           | singlet                 | 0.374426050641161     | singlet                 | singlet                 | -0.65760195       | singlet                         |  465_466                             |
    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+
    | AAACCTGAGTAATCCC-1      | singlet                 | 660_661                         | 4                       | 660_661                           | singlet                 | 0.247842972104563     | singlet                 | singlet                 | -3.5948637        | singlet                         |  660_661                             |
    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+
    | AAACCTGAGTAGCCGA-1      | doublet                 | doublet                         | unassigned              | doublet                           | doublet                 | 0.342998285281922     | singlet                 | singlet                 | -0.50507957       | doublet                         |  doublet                             |
    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+
    | ...                     | ...                     | ...                             | ...                     | ...                               | ...                     | ...                   | ...                     | ...                     | ...               | ...                             | ...                                  |
    +-------------------------+-------------------------+---------------------------------+-------------------------+-----------------------------------+-------------------------+-----------------------+-------------------------+-------------------------+-------------------+---------------------------------+--------------------------------------+

