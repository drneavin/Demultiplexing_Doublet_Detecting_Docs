.. _Scrublet-docs:


Scrublet
===========================

.. _Scrublet: https://github.com/swolock/scrublet
.. _publication: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03224-8

Scrublet_ is a transcription-based doublet detecting software.
We have provided a wrapper script that enables Scrublet_ to be easily run from the command line but we also provide example code so that users can run manually as well depending on their data.



Data
----
This is the data that you will need to have prepare to run Scrublet_:

.. admonition:: Required
  :class: important

  - A counts matrix (``$COUNTS``)
  
    - Scrublet_ expects counts to be in the cellranger output format either as

      - h5 file (``filtered_feature_bc_matrix.h5``) 
      
        **or** 
      
      - matrix directory (directory containing ``barcodes.tsv``, ``genes.tsv`` and ``matrix.mtx`` **or** ``barcodes.tsv.gz``, ``features.tsv.gz`` and ``matrix.mtx.gz``)

      - If you don't have your data in this format, you can run Scrublet_ manually in python and load the data in using a method of your choosing.

.. admonition:: Optional

  - Output directory (``$SCRUBLET_OUTDIR``)

    - If you don't provide an ``$SCRUBLET_OUTDIR``, the results will be written to the present working directory.

  - Filtered barcode file

    - A list of barcodes that are a subset of the barcodes in your h5 or matrix.mtx files. This is useful if you have run other QC softwares such as `CellBender <https://cellbender.readthedocs.io/en/stable/index.html>`__ or `DropletQC <https://github.com/powellgenomicslab/DropletQC>`__ to remove empty droplets or droplets with damaged cells.

    - Expectation is that there is no header in this file


Run Scrublet
---------------------
.. admonition:: :octicon:`stopwatch` Expected Resource Usage
  :class: note

  ~1min using a total of 15Mb memory when using 2 thread for the full :ref:`Test Dataset <TestData>` which contains ~20,982 droplets of 13 multiplexed donors,

You can either run Scrublet_ with the wrapper script we have provided or you can run it manually if you would prefer to alter more parameters.
In addition, we provide an example for filtering the single cell matrix to a subsetted list of barcodes 
 
.. admonition:: Note

  It is a good idea to try multiple different percentile variable numbers. We typically try, 80, 85, 90 and 95. 
  Then we choose the one that has the best defined bimodal distribution based on the ``doublet_score_histogram.png`` (see :ref:`Scrublet Results and Interpretation <scrublet-results>` for details).

.. tabs::

  .. tab:: With Wrapper Script

    First, let's assign the variables that will be used to execute each step.

    .. admonition:: Example Variable Settings
      :class: grey

      Below is an example of the variables that we can set up to be used in the command below.
      These are files provided as a :ref:`test dataset <TestData>` available in the :ref:`Data Preparation Documentation <DataPrep-docs>`
      Please replace paths with the full path to data on your system.

      .. code-block:: bash

        SCRUBLET_OUTDIR=/path/to/output/scrublet
        COUNTS=/path/to/TestData4PipelineFull/test_dataset/outs/filtered_gene_bc_matrices/Homo_sapiens_GRCh38p10/


    To run Scrublet_ with our wrapper script, simply execute the following in your shell:

    .. code-block:: bash

      singularity exec Demuxafy.sif Scrublet.py -m $COUNTS -o $SCRUBLET_OUTDIR

    .. admonition:: HELP! It says my file/directory doesn't exist!
      :class: dropdown

      If you receive an error indicating that a file or directory doesn't exist but you are sure that it does, this is likely an issue arising from Singularity.
      This is easy to fix.
      The issue and solution are explained in detail in the :ref:`Notes About Singularity Images <Singularity-docs>`
      
    To see all the parameters that this wrapper script will accept, run:

    .. code-block:: bash

      singularity exec Demuxafy.sif Scrublet.py -h



      usage: Scrublet.py [-h] -m COUNTS_MATRIX [-b BARCODES] [-f FILTERED_BARCODES]
                   [-r SIM_DOUBLET_RATIO] [-c MIN_COUNTS] [-e MIN_CELLS]
                   [-v MIN_GENE_VARIABILITY_PCTL] [-p N_PRIN_COMPS]
                   [-t SCRUBLET_DOUBLET_THRESHOLD] [-o OUTDIR]

      wrapper for scrublet for doublet detection of transcriptomic data.

      optional arguments:
        -h, --help            show this help message and exit
        -m COUNTS_MATRIX, --counts_matrix COUNTS_MATRIX
                              cell ranger counts matrix directory containing matrix files or full path to matrix.mtx. Can also also provide the 10x h5.
        -b BARCODES, --barcodes BARCODES
                              barcodes.tsv or barcodes.tsv.gz from cellranger
        -f FILTERED_BARCODES, --filtered_barcodes FILTERED_BARCODES
                              File containing a filtered list of droplet barcodes.
                              This may be used if you want to use a filtered list of
                              barcodes for doublet detection (ie need to remove
                              droplets that are empty or high in ambient RNA).
        -r SIM_DOUBLET_RATIO, --sim_doublet_ratio SIM_DOUBLET_RATIO
                              Number of doublets to simulate relative to the number of observed transcriptomes.
        -c MIN_COUNTS, --min_counts MIN_COUNTS
                              Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than min_cells are excluded.
        -e MIN_CELLS, --min_cells MIN_CELLS
                              Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than are excluded.
        -v MIN_GENE_VARIABILITY_PCTL, --min_gene_variability_pctl MIN_GENE_VARIABILITY_PCTL
                              Used for gene filtering prior to PCA. Keep the most highly variable genes in the top min_gene_variability_pctl percentile), as measured by the v-statistic [Klein et al., Cell 2015].
        -p N_PRIN_COMPS, --n_prin_comps N_PRIN_COMPS
                              Number of principal components used to embed the transcriptomes priorto k-nearest-neighbor graph construction.
        -t SCRUBLET_DOUBLET_THRESHOLD, --scrublet_doublet_threshold SCRUBLET_DOUBLET_THRESHOLD
                              Manually Set the scrublet doublet threshold location. For running a second time if scrublet incorrectly places the threshold the first time
        -o OUTDIR, --outdir OUTDIR
                              The output directory


  .. tab:: Run in python

    To run Scrublet_ manually, first start python from the singularity image (all the required software have been provided in the image)

    .. code-block:: bash

      singularity exec Demuxafy.sif python

  
    Now, python will open in your terminal and you can run the Scrublet_ code. 
    Here is an example:

    .. code-block:: python

      import sys
      import os
      import scrublet as scr
      import scipy.io
      import matplotlib
      matplotlib.use('AGG')
      import matplotlib.pyplot as plt
      import numpy as np
      import pandas as pd
      import umap
      import numba
      import numba.typed

      # Get path of mods directory from current script directory
      mods_path = "/opt/Demultiplexing_Doublet_Detecting_Docs/mods" ## Do not change - this is the path to the mods folder in the singularity image with custom script for loading 10x data in python
      sys.path.append(mods_path)
      import read10x

      ## Set up parameters and variables ##
      counts_matrix_dir = "/path/to/counts/matrix/dir/" ## Change this based on the path on your system
      outdir = "/path/to/doublet/detection/outdir" ## Change this based on the path on your system

      if not os.path.isdir(outdir):
        os.mkdir(outdir)


      plt.rc('font', size=14)
      plt.rcParams['pdf.fonttype'] = 42

      ## Basic run with scrublet
      counts_matrix = read10x.import_cellranger_mtx(counts_matrix_dir) ## or scanpy.read_10x_h5(counts_matrix_dir)

      try:
        barcodes_df = read10x.read_barcodes(counts_matrix_dir + "/barcodes.tsv.gz")
      except:
        try:
          barcodes_df = read10x.read_barcodes(counts_matrix_dir + "/barcodes.tsv")
        except:
          print("No barcode file in provided counts matrix directory. Please double check the directory or provide the full path to the barcode file to use.")



      dbl_rate = counts_matrix.shape[0]/1000 * 0.008 ## This is the calculation for 10x doublet rate but will be different for other platforms
      print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
      scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=dbl_rate, sim_doublet_ratio = 2)
      doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3, 
                                                                min_cells=3, 
                                                                min_gene_variability_pctl=85, 
                                                                n_prin_comps=30)


      ### Plotting and saving
      scrub.plot_histogram();
      plt.savefig(os.path.join(outdir,'doublet_score_histogram.png'))
      print('Running UMAP...')
      scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
      print('Done.')
      scrub.plot_embedding('UMAP', order_points=True);
      plt.savefig(os.path.join(outdir,'UMAP.png'))

      results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
      scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
      dataframe = pd.concat([barcodes_df, results, scores], axis=1)
      dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
      dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")

      dataframe.to_csv(os.path.join(outdir,'scrublet_results.tsv'), sep = "\t", index = False)


      ### Make summary of singlets and doublets and write to file ###
      summary = pd.DataFrame(dataframe.scrublet_DropletType.value_counts())
      summary.index.name = 'Classification'
      summary.reset_index(inplace=True)
      summary = summary.rename({'scrublet_DropletType': 'Droplet N'}, axis=1)

      summary.to_csv(os.path.join(outdir,'scrublet_summary.tsv'), sep = "\t", index = False)



  .. tab:: Run in python with filtered barcodes

    To run Scrublet_ manually, first start python from the singularity image (all the required software have been provided in the image)

    .. code-block:: bash

      singularity exec Demuxafy.sif python

  
    Now, python will open in your terminal and you can run the Scrublet_ code. 
    Here is an example:

    .. code-block:: python

      import sys
      import os
      import scrublet as scr
      import scipy.io
      import matplotlib
      matplotlib.use('AGG')
      import matplotlib.pyplot as plt
      import numpy as np
      import pandas as pd
      import umap
      import numba
      import numba.typed

      # Get path of mods directory from current script directory
      mods_path = "/opt/Demultiplexing_Doublet_Detecting_Docs/mods" ## Do not change - this is the path to the mods folder in the singularity image with custom script for loading 10x data in python
      sys.path.append(mods_path)
      import read10x

      ## Set up parameters and variables ##
      counts_matrix_dir = "/path/to/counts/matrix/dir/" ## Change this based on the path on your system
      outdir = "/path/to/doublet/detection/outdir" ## Change this based on the path on your system
      filtered_barcodes = "/path/to/filtered/barcodes/file.tsv" ## Change this based on the path on your system

      if not os.path.isdir(outdir):
        os.mkdir(outdir)


      plt.rc('font', size=14)
      plt.rcParams['pdf.fonttype'] = 42

      ## Basic run with scrublet
      counts_matrix = read10x.import_cellranger_mtx(counts_matrix_dir) ## or scanpy.read_10x_h5(counts_matrix_dir)

      try:
        barcodes_df = read10x.read_barcodes(counts_matrix_dir + "/barcodes.tsv.gz")
      except:
        try:
          barcodes_df = read10x.read_barcodes(counts_matrix_dir + "/barcodes.tsv")
        except:
          print("No barcode file in provided counts matrix directory. Please double check the directory or provide the full path to the barcode file to use.")

      ## Read in the barcodes to filter by and filter the matrix
      barcodes_filtered_df = read10x.read_barcodes(args.filtered_barcodes)

      counts_matrix = counts_matrix[barcodes_df['Barcode'].isin(barcodes_filtered_df['Barcode'])]


      dbl_rate = counts_matrix.shape[0]/1000 * 0.008 ## This is the calculation for 10x doublet rate but will be different for other platforms
      print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
      scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=dbl_rate, sim_doublet_ratio = 2)
      doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3, 
                                                                min_cells=3, 
                                                                min_gene_variability_pctl=85, 
                                                                n_prin_comps=30)


      ### Plotting and saving
      scrub.plot_histogram();
      plt.savefig(os.path.join(outdir,'doublet_score_histogram.png'))
      print('Running UMAP...')
      scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
      print('Done.')
      scrub.plot_embedding('UMAP', order_points=True);
      plt.savefig(os.path.join(outdir,'UMAP.png'))

      results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
      scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
      dataframe = pd.concat([barcodes_df, results, scores], axis=1)
      dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
      dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")

      dataframe.to_csv(os.path.join(outdir,'scrublet_results.tsv'), sep = "\t", index = False)


      ### Make summary of singlets and doublets and write to file ###
      summary = pd.DataFrame(dataframe.scrublet_DropletType.value_counts())
      summary.index.name = 'Classification'
      summary.reset_index(inplace=True)
      summary = summary.rename({'scrublet_DropletType': 'Droplet N'}, axis=1)

      summary.to_csv(os.path.join(outdir,'scrublet_summary.tsv'), sep = "\t", index = False)

  .. admonition:: HELP! I'm getting an error about 'rdist'?
    :class: dropdown

    If you receive an error similar to the following:

    ```
    RuntimeError: cannot cache function 'rdist': no locator available for file '/opt/conda/envs/py37/lib/python3.7/site-packages/umap/layouts.py'
    ```

    you may have to set the ``NUMBA_CACHE_DIR`` variable


.. _scrublet-results:

Scrublet Results and Interpretation
-------------------------------------------
After running the Scrublet_, you will have four files in the ``$SCRUBLET_OUTDIR``:

.. code-block::

  /path/to/output/scrublet
  ├── doublet_score_histogram.png
  ├── scrublet_results.tsv
  ├── scrublet_summary.tsv
  └── UMAP.png

We have found these to be the most helpful:

- ``scrublet_summary.tsv``

  - A summary of the number of singlets and doublets predicted by Scrublet_.

  +------------------------------+-----------+
  | scrublet_DropletType         | Droplet N |
  +==============================+===========+
  | doublet                      | 1851      |
  +------------------------------+-----------+
  | singlet                      | 19131     |
  +------------------------------+-----------+

    - To check whether the number of doublets identified by Scrublet_ is consistent with the expected doublet rate expected based on the number of droplets that you captured, you can use our `Expected Doublet Estimation Calculator <test.html>`__.

- ``scrublet_results.tsv``

    +------------------------+-----------------------------+-----------------+
    | Barcode                | scrublet_DropletType        | scrublet_Scores |
    +========================+=============================+=================+
    | AAACCTGAGATAGCAT-1     | singlet                     | 0.0545          |
    +------------------------+-----------------------------+-----------------+
    | AAACCTGAGCAGCGTA-1     | singlet                     | 0.1179          |
    +------------------------+-----------------------------+-----------------+
    | AAACCTGAGCGATGAC-1     | singlet                     | 0.1356          |
    +------------------------+-----------------------------+-----------------+
    | AAACCTGAGCGTAGTG-1     | singlet                     | 0.0844          |
    +------------------------+-----------------------------+-----------------+
    | AAACCTGAGGAGTTTA-1     | singlet                     | 0.0958          |
    +------------------------+-----------------------------+-----------------+
    | AAACCTGAGGCTCATT-1     | singlet                     | 0.1329          |
    +------------------------+-----------------------------+-----------------+
    | AAACCTGAGGGCACTA-1     | doublet                     | 0.4474          |
    +------------------------+-----------------------------+-----------------+
    | ...                    | ...                         | ...             |
    +------------------------+-----------------------------+-----------------+

- ``doublet_score_histogram.png``

  - This is the method that Scrublet_ uses to identify doublets - it assumes a bimodal distribution of doublet scores. Those droplets with lower scores should be singlets and those with higher scores should be doublets. It identifies the correct threshold by identifying the minimum of the bimodal distribution of simulated doublets (right).

  - However, sometimes there is not a good bimodal distribution and sometimes you will have to set the threshold manually.

  - Here is an example of a good distribution (left) and a bad distribution (left)

    +--------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
    | Good Distribution                                                                                                  | Bad Distribution                                                                                                 |
    +====================================================================================================================+==================================================================================================================+
    | .. figure:: https://user-images.githubusercontent.com/44268007/104436850-016db600-55de-11eb-8f75-229338f7bac7.png  | .. figure:: https://user-images.githubusercontent.com/44268007/88889203-ed780700-d27e-11ea-9104-60d7015f2510.png |
    +--------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+

    - In the case of the left sample, we would rerun with different parameters to try to get a better distribution and possibly manually set the threshold to ~0.2 depending on the results. In the event that we can't achieve a clear bimodal distribution, we don't use scrublet for doublet detecting.

Merging Results with Other Software Results
--------------------------------------------
We have provided a script that will help merge and summarize the results from multiple softwares together.
See :ref:`Combine Results <Combine-docs>`.


Citation
--------
If you used the Demuxafy platform for analysis, please reference our publication_ as well as `Scrublet <https://www.cell.com/cell-systems/pdfExtended/S2405-4712(18)30474-5>`__.