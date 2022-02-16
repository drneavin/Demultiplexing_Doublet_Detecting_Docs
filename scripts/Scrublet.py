#!/usr/bin/env python
import argparse
import sys
import os

parser = argparse.ArgumentParser(
    description="wrapper for scrublet for doublet detection of transcriptomic data.")
parser.add_argument("-m", "--counts_matrix", required = True, help = "cell ranger counts matrix directory containing matrix files or full path to matrix.mtx. Can also also provide the 10x h5.")
parser.add_argument("-b", "--barcodes", required = False, help = "barcodes.tsv or barcodes.tsv.gz from cellranger")
parser.add_argument("-r", "--sim_doublet_ratio", required = False, default = 2, type = int, help = "Number of doublets to simulate relative to the number of observed transcriptomes.")
parser.add_argument("-c", "--min_counts", required = False, default = 3, type = int, help = "Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than min_cells are excluded.")
parser.add_argument("-e", "--min_cells", required = False, default = 3, type = int, help = "Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than are excluded.")
parser.add_argument("-v", "--min_gene_variability_pctl", required = False, default = 85, type = int, help = "Used for gene filtering prior to PCA. Keep the most highly variable genes in the top min_gene_variability_pctl percentile), as measured by the v-statistic [Klein et al., Cell 2015].")
parser.add_argument("-p", "--n_prin_comps", required = False, default = 30, type = int, help = "Number of principal components used to embed the transcriptomes priorto k-nearest-neighbor graph construction.")
parser.add_argument("-t", "--scrublet_doublet_threshold", required = False, default = None, type = float, help = "Manually Set the scrublet doublet threshold location. For running a second time if scrublet incorrectly places the threshold the first time")
parser.add_argument("-o", "--outdir", required = False, default = os.getcwd(), help = "The output directory")
args = parser.parse_args()

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

if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)


plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

## Basic run with scrublet
if os.path.exists(args.counts_matrix):
    if args.counts_matrix.endswith(".h5"):
        counts_matrix = scanpy.read_10x_h5(args.counts_matrix)
    else:
        counts_matrix = read10x.import_cellranger_mtx(args.counts_matrix)
else:
    print("Couldn't find the counts file " + args.counts_matrix)


if args.barcodes is None:
	if os.path.exists(os.path.join(args.counts_matrix, "barcodes.tsv.gz")):
		print("Reading in barcodes file")
		barcodes_df = read10x.read_barcodes(os.path.join(args.counts_matrix ,"barcodes.tsv.gz"))
	elif os.path.exists(os.path.join(args.counts_matrix, "barcodes.tsv")):
		print("Reading in barcodes file")
		barcodes_df = read10x.read_barcodes(os.path.join(args.counts_matrix ,"barcodes.tsv"))
	else:
		print("No barcode file in provided counts matrix directory")
else:
	print("Reading in barcodes file")
	barcodes_df = read10x.read_barcodes(args.barcodes)



dbl_rate = counts_matrix.shape[0]/1000 * 0.008
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=dbl_rate, sim_doublet_ratio = args.sim_doublet_ratio)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=args.min_counts, 
                                                          min_cells=args.min_cells, 
                                                          min_gene_variability_pctl=args.min_gene_variability_pctl, 
                                                          n_prin_comps=args.n_prin_comps)


if args.scrublet_doublet_threshold is None:
  ### Plotting and saving
  scrub.plot_histogram();
  plt.savefig(os.path.join(args.outdir,'doublet_score_histogram.png'))
  print('Running UMAP...')
  scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
  print('Done.')
  scrub.plot_embedding('UMAP', order_points=True);
  plt.savefig(os.path.join(args.outdir,'UMAP.png'))

  results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
  scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
  dataframe = pd.concat([barcodes_df, results, scores], axis=1)
  dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
  dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")

  dataframe.to_csv(os.path.join(args.outdir,'scrublet_results.tsv'), sep = "\t", index = False)

else:
  print(args.scrublet_doublet_threshold)
  scrub.call_doublets(threshold=args.scrublet_doublet_threshold)
  ### Plotting and saving
  scrub.plot_histogram(); 
  plt.savefig(os.path.join(args.outdir,'doublet_score_histogram_manual_threshold.png'))
  print('Running UMAP...')
  scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
  print('Done.')
  scrub.plot_embedding('UMAP', order_points=True);
  plt.savefig(os.path.join(args.outdir,'UMAP_manual_threshold.png'))

  results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
  scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
  dataframe = pd.concat([barcodes_df, results, scores], axis=1)
  dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
  dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")
  
  dataframe.to_csv(os.path.join(args.outdir,'scrublet_doublets_singlets.tsv'), sep = "\t", index = False)


### Make summary of singlets and doublets and write to file ###
summary = pd.DataFrame(dataframe.scrublet_DropletType.value_counts())
summary.index.name = 'Classification'
summary.reset_index(inplace=True)
summary = summary.rename({'scrublet_DropletType': 'Droplet N'}, axis=1)

summary.to_csv(os.path.join(args.outdir,'scrublet_summary.tsv'), sep = "\t", index = False)

