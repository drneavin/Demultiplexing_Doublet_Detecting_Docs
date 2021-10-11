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




Citation
--------
If you used this workflow for analysis, please reference our paper (REFERENCE).