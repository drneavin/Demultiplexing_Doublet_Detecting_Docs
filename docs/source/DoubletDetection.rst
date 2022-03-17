.. _DoubletDetection-docs:

DoubletDetection
===========================

.. _DoubletDetection: https://github.com/JonathanShor/DoubletDetection
.. _preprint: https://www.biorxiv.org/content/10.1101/2022.03.07.483367v1


DoubletDetection_ is a transcription-based doublet detection software.
This was one of the better-performing doublet detecting softwares that we identified in our paper (**CITE**) and it is also relatively fast to run.
We have provided a wrapper script that enables DoubletDetection_ to be easily run from the command line but we also provide example code so that users can run manually as well depending on their data.



Data
----
This is the data that you will need to have prepare to run DoubletDetection_:

.. admonition:: Required
  :class: important

  - A counts matrix (``$COUNTS``)
  
    - DoubletDetection_ expects counts to be in the cellranger output format either as
    
      - h5 file (``filtered_feature_bc_matrix.h5``) 
      
        **or** 
      
      - matrix directory (directory containing ``barcodes.tsv``, ``genes.tsv`` and ``matrix.mtx`` **or** ``barcodes.tsv.gz``, ``features.tsv.gz`` and ``matrix.mtx.gz``)

      - If you don't have your data in either of these formats, you can run DoubletDetection_ manually in python and load the data in using a method of your choosing.

.. admonition:: Optional

  - Output directory (``$DOUBLETDETECTION_OUTDIR``)

    - If you don't provide an ``$DOUBLETDETECTION_OUTDIR``, the results will be written to the present working directory.


Run DoubletDetection
---------------------
You can either run DoubletDetection_ with the wrapper script we have provided or you can run it manually if you would prefer to alter more parameters.

.. tabs::

  .. tab:: With Wrapper Script

    .. code-block:: bash

      singularity exec Demuxafy.sif DoubletDetection.py -m $COUNTS -o $DOUBLETDETECTION_OUTDIR

    To see all the parameters that this wrapper script will accept, run:

    .. code-block:: bash

      singularity exec Demuxafy.sif DoubletDetection.py -h

      usage: DoubletDetection.py [-h] -m COUNTS_MATRIX [-b BARCODES] [-o OUTDIR] [-i N_ITERATIONS] [-p PHENOGRAPH] [-s STANDARD_SCALING] [-t P_THRESH] [-v VOTER_THRESH]

      wrapper for DoubletDetection for doublet detection from transcriptomic data.

      optional arguments:
        -h, --help            show this help message and exit
        -m COUNTS_MATRIX, --counts_matrix COUNTS_MATRIX
                              cell ranger counts matrix directory containing matrix files or full path to matrix.mtx. Can also also provide the 10x h5.
        -b BARCODES, --barcodes BARCODES
                              File containing droplet barcodes. Use barcodes from provided 10x dir by default.
        -o OUTDIR, --outdir OUTDIR
                              The output directory; default is current working directory
        -i N_ITERATIONS, --n_iterations N_ITERATIONS
                              Number of iterations to use; default is 50
        -p PHENOGRAPH, --phenograph PHENOGRAPH
                              Whether to use phenograph (True) or not (False); default is False
        -s STANDARD_SCALING, --standard_scaling STANDARD_SCALING
                              Whether to use standard scaling of normalized count matrix prior to PCA (True) or not (False); default is True
        -t P_THRESH, --p_thresh P_THRESH
                              P-value threshold for doublet calling; default is 1e-16
        -v VOTER_THRESH, --voter_thresh VOTER_THRESH
                              Voter threshold for doublet calling; default is 0.5

  .. tab:: Run in python

    To run DoubletDetection_ manually, first start python from the singularity image (all the required software have been provided in the image)

    .. code-block:: bash

      singularity exec Demuxafy.sif python

    Now, python will open in your terminal and you can run the DoubletDetection_ code. 
    Here is an example:

    .. code-block:: python

      import os
      import numpy as np
      import doubletdetection
      import tarfile
      import matplotlib
      matplotlib.use('PDF')
      import matplotlib.pyplot as plt
      import sys
      import pandas as pd

      # Load read10x function from mods directory

      mods_path = "/opt/Demultiplexing_Doublet_Detecting_Docs/mods" ## custom script for loading 10x data in python
      sys.path.append(mods_path)
      import read10x

      ### Set up parameters and variables ###
      counts_matrix = "/path/to/counts/matrix.mtx"
      outdir = "/path/to/doublet/detection/outdir"


      if not os.path.isdir(outdir):
      	os.mkdir(outdir)


      ### Read in data ###
      raw_counts = read10x.import_cellranger_mtx(counts_matrix)

      try:
        barcodes_df = read10x.read_barcodes(counts_matrix + "/barcodes.tsv.gz")
      except:
        try:
          barcodes_df = read10x.read_barcodes(counts_matrix + "/barcodes.tsv")
        except:
          print("No barcode file in provided counts matrix directory. Please double check the directory or provide the full path to the barcode file to use.")

      print('Counts matrix shape: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

      # Remove columns with all 0s
      zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
      raw_counts = raw_counts[:, ~zero_genes]
      print('Counts matrix shape after removing unexpressed genes: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

      clf = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose = True)
      doublets = clf.fit(raw_counts).predict(p_thresh=1e-16, voter_thresh=50)

      results = pd.Series(doublets, name="DoubletDetection_DropletType")
      dataframe = pd.concat([barcodes_df, results], axis=1)
      dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(1.0, "doublet")
      dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(0.0, "singlet")

      dataframe.to_csv(os.path.join(outdir,'DoubletDetection_doublets_singlets.tsv'), sep = "\t", index = False)


      ### Figures ###
      doubletdetection.plot.convergence(clf, save=os.path.join(outdir,'convergence_test.pdf'), show=False, p_thresh=1e-16, voter_thresh=0.5)

      f3 = doubletdetection.plot.threshold(clf, save=os.path.join(outdir,'threshold_test.pdf'), show=False, p_step=6)


      ### Make summary of singlets and doublets and write to file ###
      summary = pd.DataFrame(dataframe.DoubletDetection_DropletType.value_counts())
      summary.index.name = 'Classification'
      summary.reset_index(inplace=True)
      summary = summary.rename({'DoubletDetection_DropletType': 'Droplet N'}, axis=1)

      summary.to_csv(os.path.join(outdir,'DoubletDetection_summary.tsv'), sep = "\t", index = False)




DoubletDetection Results and Interpretation
-------------------------------------------
After running the DoubletDetection_, you will have multiple files in the ``$DOUBLETDETECTION_OUTDIR``:

.. code-block:: bash

  .
  ├── convergence_test.pdf
  ├── DoubletDetection_doublets_singlets.tsv
  ├── DoubletDetection_summary.tsv
  └── threshold_test.pdf

We have found these to be the most helpful:

- ``DoubletDetection_summary.tsv``

  - A summary of the number of singlets and doublets predicted by DoubletDetection_.

  +------------------------------+-----------+
  | DoubletDetection_DropletType | Droplet N |
  +==============================+===========+
  | doublet                      | 2594      |
  +------------------------------+-----------+
  | singlet                      | 18388     |
  +------------------------------+-----------+

    - To check whether the number of doublets identified by DoubletDetection_ is consistent with the expected doublet rate expected based on the number of droplets that you captured, you can use our `Expected Doublet Estimation Calculator <test.html>`__.

- ``DoubletDetection_doublets_singlets.tsv``

  - The per-barcode singlet and doublet classification from DoubletDetection_.

    +------------------------+-----------------------------+
    | Barcode                | DoubletDetection_DropletType|
    +========================+=============================+
    | AAACCTGAGATAGCAT-1     | singlet                     |
    +------------------------+-----------------------------+
    | AAACCTGAGCAGCGTA-1     | singlet                     |
    +------------------------+-----------------------------+
    | AAACCTGAGCGATGAC-1     | singlet                     |
    +------------------------+-----------------------------+
    | AAACCTGAGCGTAGTG-1     | singlet                     |
    +------------------------+-----------------------------+
    | AAACCTGAGGAGTTTA-1     | singlet                     |
    +------------------------+-----------------------------+
    | AAACCTGAGGCTCATT-1     | singlet                     |
    +------------------------+-----------------------------+
    | AAACCTGAGGGCACTA-1     | singlet                     |
    +------------------------+-----------------------------+
    | ...                    | ...                         |
    +------------------------+-----------------------------+

- ``convergence_test.pdf``

  - The expectation is that after multiple rounds, the expected number of doublets will converge. If that is not the case, we suggest that you run DoubletDetection for more iterations (try 150, or even 250 if that isn't convincing).

  - Here are two figures - one of a sample that came to convergence after 50 iterations (left) and one that did not (right)

    +--------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
    | Good Converged                                                                                                     | Bad Convergence                                                                                                  |
    +====================================================================================================================+==================================================================================================================+
    | .. figure:: https://user-images.githubusercontent.com/44268007/104434976-ccf8fa80-55db-11eb-9f30-00f71e4592d4.png  | .. figure:: https://user-images.githubusercontent.com/44268007/95423527-f545dd00-098c-11eb-8a48-1ca6bb507151.png |
    +--------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+


Merging Results with Other Software Results
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.

Citation
--------
If you used the Demuxafy platform for analysis, please reference our preprint_ as well as `DoubletDetection <https://zenodo.org/record/4359992>`__.