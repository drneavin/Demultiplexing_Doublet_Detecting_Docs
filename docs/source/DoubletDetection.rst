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

  - Filtered barcode file

    - A list of barcodes that are a subset of the barcodes in your h5 or matrix.mtx files. This is useful if you have run other QC softwares such as `CellBender <https://cellbender.readthedocs.io/en/stable/index.html>`__ or `DropletQC <https://github.com/powellgenomicslab/DropletQC>`__ to remove empty droplets or droplets with damaged cells.

    - Expectation is that there is no header in this file



Run DoubletDetection
---------------------
.. admonition:: |:stopwatch:| Expected Resource Usage
  :class: note

  ~1h using a total of 15Gb memory when using 2 thread for the full :ref:`Test Dataset <TestData>` which contains ~20,982 droplets of 13 multiplexed donors,

You can either run DoubletDetection_ with the wrapper script we have provided or you can run it manually if you would prefer to alter more parameters.
In addition, we provide an example for filtering the single cell matrix to a subsetted list of barcodes 

.. tabs::

  .. tab:: With Wrapper Script

    First, let's assign the variables that will be used to execute each step.

    .. admonition:: Example Variable Settings
      :class: grey

      Below is an example of the variables that we can set up to be used in the command below.
      These are files provided as a :ref:`test dataset <TestData>` available in the :ref:`Data Preparation Documentation <DataPrep-docs>`
      Please replace paths with the full path to data on your system.

      .. code-block:: bash

        DOUBLETDETECTION_OUTDIR=/path/to/output/DoubletDetection
        COUNTS=/path/to/TestData4PipelineFull/test_dataset/outs/filtered_gene_bc_matrices/Homo_sapiens_GRCh38p10/


    .. code-block:: bash

      singularity exec Demuxafy.sif DoubletDetection.py -m $COUNTS -o $DOUBLETDETECTION_OUTDIR

    .. admonition:: HELP! It says my file/directory doesn't exist!
      :class: dropdown

      If you receive an error indicating that a file or directory doesn't exist but you are sure that it does, this is likely an issue arising from Singularity.
      This is easy to fix.
      The issue and solution are explained in detail in the :ref:`Notes About Singularity Images <Singularity-docs>`

    To see all the parameters that this wrapper script will accept, run:

    .. code-block:: bash

      singularity exec Demuxafy.sif DoubletDetection.py -h

      usage: DoubletDetection.py [-h] -m COUNTS_MATRIX [-b BARCODES]
                                [-f FILTERED_BARCODES] [-o OUTDIR] [-p BOOST_RATE]
                                [-c N_COMPONENTS] [-g N_TOP_VAR_GENES] [-r REPLACE]
                                [-a CLUSTERING_ALGORITHM] [-k CLUSTERING_KWARGS]
                                [-i N_ITERATIONS] [-e PSEUDOCOUNT] [-n NORMALIZER]
                                [-d RANDOM_STATE] [-s STANDARD_SCALING] [-j N_JOBS]
                                [-t P_THRESH] [-v VOTER_THRESH]

      wrapper for DoubletDetection for doublet detection from transcriptomic data.

      optional arguments:
        -h, --help            show this help message and exit
        -m COUNTS_MATRIX, --counts_matrix COUNTS_MATRIX
                              cell ranger counts matrix directory containing matrix
                              files or full path to matrix.mtx. Can also also
                              provide the 10x h5.
        -b BARCODES, --barcodes BARCODES
                              File containing droplet barcodes. Use barcodes from
                              provided 10x dir by default.
        -f FILTERED_BARCODES, --filtered_barcodes FILTERED_BARCODES
                              File containing a filtered list of droplet barcodes.
                              This may be used if you want to use a filtered list of
                              barcodes for doublet detection (ie need to remove
                              droplets that are empty or high in ambient RNA).
        -o OUTDIR, --outdir OUTDIR
                              The output directory; default is current working
                              directory
        -p BOOST_RATE, --boost_rate BOOST_RATE
                              Proportion of cells used to generate synthetic
                              doublets; default is 0.25.
        -c N_COMPONENTS, --n_components N_COMPONENTS
                              Number of principal components to use; default is 30.
        -g N_TOP_VAR_GENES, --n_top_var_genes N_TOP_VAR_GENES
                              Number of top variable genes to use; default is 1000.
        -r REPLACE, --replace REPLACE
                              Whether to replace cells when generating synthetic
                              doublets; default is False.
        -a CLUSTERING_ALGORITHM, --clustering_algorithm CLUSTERING_ALGORITHM
                              Which clustering algorithm to use; default is
                              'phenograph'
        -k CLUSTERING_KWARGS, --clustering_kwargs CLUSTERING_KWARGS
                              Keyword arguments to pass to clustering algorithm;
                              default is None.
        -i N_ITERATIONS, --n_iterations N_ITERATIONS
                              Number of iterations to use; default is 50
        -e PSEUDOCOUNT, --pseudocount PSEUDOCOUNT
                              Pseudocount used to normalize counts; default is 0.1.
        -n NORMALIZER, --normalizer NORMALIZER
                              Method for raw counts normalization; default is None.
        -d RANDOM_STATE, --random_state RANDOM_STATE
                              Number to use to seed random state for PCA; default is
                              0.
        -s STANDARD_SCALING, --standard_scaling STANDARD_SCALING
                              Whether to use standard scaling of normalized count
                              matrix prior to PCA (True) or not (False); default is
                              True
        -j N_JOBS, --n_jobs N_JOBS
                              Number of jobs to to use; default is 1
        -t P_THRESH, --p_thresh P_THRESH
                              P-value threshold for doublet calling; default is
                              1e-16
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

      mods_path = "/opt/Demultiplexing_Doublet_Detecting_Docs/mods" ## Do not change - this is the path to the mods folder in the singularity image with custom script for loading 10x data in python
      sys.path.append(mods_path)
      import read10x

      ### Set up parameters and variables ###
      counts_matrix = "/path/to/counts/matrix.mtx" ## Change this based on the path on your system
      outdir = "/path/to/doublet/detection/outdir" ## Change this based on the path on your system


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

      clf = doubletdetection.BoostClassifier(n_iters=50, clustering_algorithm='phenograph', standard_scaling=True, verbose = True)
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




  .. tab:: Run in python with filtered barcodes

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

      mods_path = "/opt/Demultiplexing_Doublet_Detecting_Docs/mods" ## Do not change - this is the path to the mods folder in the singularity image with custom script for loading 10x data in python
      sys.path.append(mods_path)
      import read10x

      ### Set up parameters and variables ###
      counts_matrix = "/path/to/counts/matrix.mtx" ## Change this based on the path on your system
      outdir = "/path/to/doublet/detection/outdir" ## Change this based on the path on your system
      filtered_barcodes = "/path/to/filtered/barcodes/file.tsv" ## Change this based on the path on your system


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


      ## Read in the barcodes to filter by and filter the matrix
      barcodes_filtered_df = read10x.read_barcodes(args.filtered_barcodes)

      raw_counts = raw_counts[barcodes_df['Barcode'].isin(barcodes_filtered_df['Barcode'])]


      clf = doubletdetection.BoostClassifier(n_iters=50, clustering_algorithm='phenograph', standard_scaling=True, verbose = True)
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
After running DoubletDetection_, you will have multiple files in the ``$DOUBLETDETECTION_OUTDIR``:

.. code-block:: bash

  /path/to/output/DoubletDetection
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