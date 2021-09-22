.. _Scrublet-docs:


Scrublet Tutorial
===========================

.. _Scrublet: https://github.com/swolock/scrublet

Scrublet_ is a transcription-based doublet detecting software.
We have provided a wrapper script that enables Scrublet_ to be easily run from the command line but we also provide example code so that users can run manually as well depending on their data.



Data
----
This is the data that you will need to have preparede to run DoubletDetection_:

.. admonition:: Required
  :class: important



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



















Citation
--------
If you used this workflow for analysis, please reference our paper (REFERENCE) as well as `Scrublet <https://www.cell.com/cell-systems/pdfExtended/S2405-4712(18)30474-5>`__.