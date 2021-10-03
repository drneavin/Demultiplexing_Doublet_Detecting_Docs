.. _Scrublet-docs:


Scrublet Tutorial
===========================

.. _Scrublet: https://github.com/swolock/scrublet

Scrublet_ is a transcription-based doublet detecting software.
We have provided a wrapper script that enables Scrublet_ to be easily run from the command line but we also provide example code so that users can run manually as well depending on their data.



Data
----
This is the data that you will need to have preparede to run Scrublet_:

.. admonition:: Required
  :class: important

  - A counts matrix (``$COUNTS``)
  
    - DoubletDetection expects counts to be in the cellranger output format (directory containint ``barcodes.tsv``, ``genes.tsv`` and ``matrix.mtx`` **or** ``barcodes.tsv.gz``, ``features.tsv.gz`` and ``matrix.mtx.gz``)

	  - If you don't have your data in this format, you can run Scrublet_ manually in python and load the data in using a method of your choosing.

.. admonition:: Optional

  - Output directory (``$OUTDIR``)

    - If you don't provide an ``$OUTDIR``, the results will be written to the present working directory.



Run Scrublet
---------------------
You can either run Scrublet_ with the wrapper script we have provided or you can run it manually if you would prefer to alter more parameters.
 
.. admonition:: Note

  It is a good idea to try multiple different percentile variable numbers. We typically try, 80, 85, 90 and 95. 
  Then we choose the one that has the best defined bimodal distribution based on the ``doublet_score_histogram.png`` (see results explanation for details).

.. tabs::

  .. tab:: With Wrapper Script

    .. code-block:: bash

      singularity exec Demuxafy.sif python scrublet.py

	  To see all the parameters that this wrapper script will accept, run:

	  ..code-block:: bash

			python scrublet.py -h


			usage: scrublet.py [-h] -m COUNTS_MATRIX [-b BARCODES] [-r SIM_DOUBLET_RATIO]
							[-c MIN_COUNTS] [-e MIN_CELLS]
							[-v MIN_GENE_VARIABILITY_PCTL] [-p N_PRIN_COMPS]
							[-t SCRUBLET_DOUBLET_THRESHOLD] [-o OUTDIR]

			wrapper for scrublet for doublet detection of transcriptomic data.

			optional arguments:
			-h, --help            show this help message and exit
			-m COUNTS_MATRIX, --counts_matrix COUNTS_MATRIX
									cell ranger counts matrix directory
			-b BARCODES, --barcodes BARCODES
									barcodes.tsv or barcodes.tsv.gz from cellranger
			-r SIM_DOUBLET_RATIO, --sim_doublet_ratio SIM_DOUBLET_RATIO
									Number of doublets to simulate relative to the number
									of observed transcriptomes.
			-c MIN_COUNTS, --min_counts MIN_COUNTS
									Used for gene filtering prior to PCA. Genes expressed
									at fewer than min_counts in fewer than min_cells are
									excluded.
			-e MIN_CELLS, --min_cells MIN_CELLS
									Used for gene filtering prior to PCA. Genes expressed
									at fewer than min_counts in fewer than are excluded.
			-v MIN_GENE_VARIABILITY_PCTL, --min_gene_variability_pctl MIN_GENE_VARIABILITY_PCTL
									Used for gene filtering prior to PCA. Keep the most
									highly variable genes in the top
									min_gene_variability_pctl percentile), as measured by
									the v-statistic [Klein et al., Cell 2015].
			-p N_PRIN_COMPS, --n_prin_comps N_PRIN_COMPS
									Number of principal components used to embed the
									transcriptomes priorto k-nearest-neighbor graph
									construction.
			-t SCRUBLET_DOUBLET_THRESHOLD, --scrublet_doublet_threshold SCRUBLET_DOUBLET_THRESHOLD
									Manually Set the scrublet doublet threshold location.
									For running a second time if scrublet incorreclty
									places the threshold the first time
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
      mods_path = "/opt/Demultiplexing_Doublet_Detecting_Docs/mods"
      sys.path.append(mods_path)
      import read10x

      ## Set up parameters and variables ##
      counts_matrix = "/path/to/counts/matrix.mtx"
      outdir = "/path/to/doublet/detection/outdir"


      os.mkdirs(outdir)


      plt.rc('font', size=14)
      plt.rcParams['pdf.fonttype'] = 42

      ## Basic run with scrublet
      counts_matrix = read10x.import_cellranger_mtx(counts_matrix)
      barcodes_df = read10x.read_barcodes(barcodes)


      dbl_rate = counts_matrix.shape[0]/1000 * 0.008
      print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
      scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=dbl_rate, sim_doublet_ratio = 2)
      doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3, 
                                                                min_cells=a3, 
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



DoubletDetection Results and Interpretation
-------------------------------------------
After running the Scrublet_, you will have multiple files in the ``$OUTDIR``.
We have found these to be the most helpful:

- ``scrublet_summary.tsv``

  - A sumamry of the number of singlets and doublets predicted by Scrublet_.

  +------------------------------+-----------+
  | scrublet_DropletType         | Droplet N |
  +==============================+===========+
  | doublet                      | 1851      |
  +------------------------------+-----------+
  | singlet                      | 19131     |
  +------------------------------+-----------+

    - To check whether the numbe of doublets identified by Scrublet_ is consistent with the expected doublet rate expected based on the number of droplets that you captured, you can use our `Expected Doublet Estimation Calculator <test.html>`__.

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




Citation
--------
If you used this workflow for analysis, please reference our paper (REFERENCE) as well as `Scrublet <https://www.cell.com/cell-systems/pdfExtended/S2405-4712(18)30474-5>`__.