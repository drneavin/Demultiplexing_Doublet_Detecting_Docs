.. _solo-docs:

Solo
===========================

.. _solo: https://github.com/calico/solo

Solo_ is a transcription-based doublet detecting software that was one of the better transcription-based doublet detecting softwares that we tested (CITATION).



Data
----
This is the data that you will need to have prepare to run Solo_:

.. admonition:: Required
  :class: important

  - Parameter json file (``$JSON``)
  
    - Solo_ has provided an :download:`example file <_download_files/solo_model.json>` that we have found to work well for most of our data.

  - Counts (``$COUNTS``)

    - This can be a h5ad file, loom file, or 10x counts matrix directory (containing ``barcodes.tsv``, ``genes.tsv`` and ``matrix.mtx`` **or** ``barcodes.tsv.gz``, ``features.tsv.gz`` and ``matrix.mtx.gz``)

  - Output directory (``$SOLO_OUTDIR``)


.. admonition:: Optional

  - Expected number of doublets (``$N_DOUB``)



Run Solo
----------------

.. code-block:: bash

  singularity exec Demuxafy.sif solo -o $SOLO_OUTDIR -e $N_DOUB -j $JSON -d $COUNTS

Solo_ also has additional parameters that can be seen with:

.. code-block:: bash

  singularity exec Demuxafy.sif solo -h 

  usage: solo [-h] -j MODEL_JSON_FILE -d DATA_PATH
            [--set-reproducible-seed REPRODUCIBLE_SEED]
            [--doublet-depth DOUBLET_DEPTH] [-g] [-a] [-o OUT_DIR]
            [-r DOUBLET_RATIO] [-s SEED] [-e EXPECTED_NUMBER_OF_DOUBLETS] [-p]
            [-recalibrate_scores] [--version]

  optional arguments:
    -h, --help            show this help message and exit
    -j MODEL_JSON_FILE    json file to pass VAE parameters (default: None)
    -d DATA_PATH          path to h5ad, loom, or 10x mtx dir cell by genes
                          counts (default: None)
    --set-reproducible-seed REPRODUCIBLE_SEED
                          Reproducible seed, give an int to set seed (default:
                          None)
    --doublet-depth DOUBLET_DEPTH
                          Depth multiplier for a doublet relative to the average
                          of its constituents (default: 2.0)
    -g                    Run on GPU (default: True)
    -a                    output modified anndata object with solo scores Only
                          works for anndata (default: False)
    -o OUT_DIR
    -r DOUBLET_RATIO      Ratio of doublets to true cells (default: 2)
    -s SEED               Path to previous solo output directory. Seed VAE
                          models with previously trained solo model. Directory
                          structure is assumed to be the same as solo output
                          directory structure. should at least have a vae.pt a
                          pickled object of vae weights and a latent.npy an
                          np.ndarray of the latents of your cells. (default:
                          None)
    -e EXPECTED_NUMBER_OF_DOUBLETS
                          Experimentally expected number of doublets (default:
                          None)
    -p                    Plot outputs for solo (default: False)
    -recalibrate_scores   Recalibrate doublet scores (not recommended anymore)
                          (default: False)
    --version             Get version of solo-sc (default: False)

If Solo_ runs correctly, you should have the following files and directory structure in your ``$SOLO_OUTDIR``:

.. code-block::

  .
  ├── classifier
  │   ├── attr.pkl
  │   ├── model_params.pt
  │   └── var_names.csv
  ├── is_doublet.csv
  ├── is_doublet.npy
  ├── is_doublet_sim.npy
  ├── latent.npy
  ├── logit_scores.csv
  ├── logit_scores.npy
  ├── logit_scores_sim.npy
  ├── no_updates_softmax_scores.csv
  ├── no_updates_softmax_scores.npy
  ├── no_updates_softmax_scores_sim.npy
  ├── preds.csv
  ├── preds.npy
  ├── smoothed_preds.npy
  ├── softmax_scores.csv
  ├── softmax_scores.npy
  └── vae
      ├── attr.pkl
      ├── model_params.pt
      └── var_names.csv


Solo Summary
^^^^^^^^^^^^^^^^
We have provided a script that will summarize the number of droplets classified as doublets and singlets by Solo_ and write it to the ``$SOLO_OUTDIR``.
This script also combines some of the Solo_ outputs into a single file that can be more easily used for downstream analyses. 
You can run this to get a fast and easy summary of your results with:

.. code-block:: bash

  singularity exec Demuxafy.sif solo_summary.py -b $BARCODES -s $SOLO_OUTDIR

If successful, you should have two new files in your ``$SOLO_OUTDIR``:

.. code-block::
  :emphasize-lines: 21,22

  .
  ├── classifier
  │   ├── attr.pkl
  │   ├── model_params.pt
  │   └── var_names.csv
  ├── is_doublet.csv
  ├── is_doublet.npy
  ├── is_doublet_sim.npy
  ├── latent.npy
  ├── logit_scores.csv
  ├── logit_scores.npy
  ├── logit_scores_sim.npy
  ├── no_updates_softmax_scores.csv
  ├── no_updates_softmax_scores.npy
  ├── no_updates_softmax_scores_sim.npy
  ├── preds.csv
  ├── preds.npy
  ├── smoothed_preds.npy
  ├── softmax_scores.csv
  ├── softmax_scores.npy
  ├── solo_results.tsv
  ├── solo_summary.tsv
  └── vae
      ├── attr.pkl
      ├── model_params.pt
      └── var_names.csv


Solo Results and Interpretation
----------------------------------------
Solo_ puts most of the results in multiple separate files. 
However, the wrapper script and the example code has some steps to combine these results together into a single file, which will likely be the most informative output.

- ``solo_summary.tsv``

  - A summary of the number of singlets and doublets predicted by Solo_.

  +-----------------+-----------+
  | Classification  | Droplet N |
  +=================+===========+
  | singlet         | 17461     |
  +-----------------+-----------+
  | doublet         | 3521      |
  +-----------------+-----------+

    - To check whether the number of doublets identified by Solo_ is consistent with the expected doublet rate expected based on the number of droplets that you captured, you can use our `Expected Doublet Estimation Calculator <test.html>`__.

- ``solo_results.tsv``

  - The per-barcode singlet and doublet classification from Solo_.

    +-------------------------+-------------------------+--------------------------+
    | Barcode                 | solo_DropletType        | solo_DropletScore        |
    +=========================+=========================+==========================+
    | AAACCTGAGATAGCAT-1      | singlet                 | -8.442187                |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGCAGCGTA-1      | singlet                 | -2.8096201               |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGCGATGAC-1      | singlet                 | -2.8949203               |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGCGTAGTG-1      | singlet                 | -5.928284                |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGGAGTTTA-1      | doublet                 | 0.2749935                |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGGCTCATT-1      | singlet                 | -5.2726507               |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGGGCACTA-1      | singlet                 | -0.65760195              |
    +-------------------------+-------------------------+--------------------------+
    | AAACCTGAGTAATCCC-1      | singlet                 | -3.5948637               |
    +-------------------------+-------------------------+--------------------------+
    | ...                     | ...                     | ...                      |
    +-------------------------+-------------------------+--------------------------+


Merging Results with Other Software Results
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.


Citation
--------
If you used the Demuxafy platform for analysis, please reference our paper (REFERENCE) as well as `solo <https://www.sciencedirect.com/science/article/pii/S2405471220301952>`__.