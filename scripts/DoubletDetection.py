#!/usr/bin/env python
import argparse
import os


parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("-m", "--counts_matrix", required = True, help = "cell ranger counts matrix directory containing matrix files or full path to matrix.mtx. Can also also provide the 10x h5.")
parser.add_argument("-b", "--barcodes", required = False, default = None, help = "File containing droplet barcodes. Use barcodes from provided 10x dir by default.")
parser.add_argument("-f", "--filtered_barcodes", required = False, default = None, help = "File containing a filtered list of droplet barcodes. This may be used if you want to use a filtered list of barcodes for doublet detection (ie need to remove droplets that are empty or high in ambient RNA).")
parser.add_argument("-o", "--outdir", required = False, default = os.getcwd(), help = "The output directory; default is current working directory")
parser.add_argument("-i", "--n_iterations", required = False, default = 50, type = int, help = "Number of iterations to use; default is 50")
parser.add_argument("-p", "--phenograph", required = False, default = False, help = "Whether to use phenograph (True) or not (False); default is False")
parser.add_argument("-s", "--standard_scaling", required = False, default = True, help = "Whether to use standard scaling of normalized count matrix prior to PCA (True) or not (False); default is True")
parser.add_argument("-t", "--p_thresh", required = False, default = 1e-16, type = float, help = "P-value threshold for doublet calling; default is 1e-16")
parser.add_argument("-v", "--voter_thresh", required = False, default = 0.5, type = float, help = "Voter threshold for doublet calling; default is 0.5")
args = parser.parse_args()

import numpy as np
import doubletdetection
import tarfile
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import scanpy

# Load read10x function from mods directory

mods_path = "/opt/Demultiplexing_Doublet_Detecting_Docs/mods"
sys.path.append(mods_path)
import read10x



if args.phenograph == 'True':
    pheno = True
elif args.phenograph == 'False':
    pheno = False
else:
    pheno = args.phenograph


if args.standard_scaling == 'True':
    standard_scaling = True
elif args.standard_scaling == 'False':
    standard_scaling = False
else:
    standard_scaling = args.standard_scaling

if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)


### Read in data ###
if os.path.exists(args.counts_matrix):
    if args.counts_matrix.endswith(".h5"):
        raw_counts = scanpy.read_10x_h5(args.counts_matrix)
    else:
        raw_counts = read10x.import_cellranger_mtx(args.counts_matrix)
else:
    print("Couldn't find the counts file " + args.counts_matrix)


if args.barcodes is None:
    if os.path.exists(os.path.join(args.counts_matrix, "barcodes.tsv.gz")):
        barcodes_df = read10x.read_barcodes(os.path.join(args.counts_matrix ,"barcodes.tsv.gz"))
    elif os.path.exists(os.path.join(args.counts_matrix, "barcodes.tsv")):
        barcodes_df = read10x.read_barcodes(os.path.join(args.counts_matrix ,"barcodes.tsv"))
    else:
        print("No barcode file in provided or couldn't find it at counts matrix directory " + args.counts_matrix)
        exit()
else:
    barcodes_df = read10x.read_barcodes(args.barcodes)

print('Counts matrix shape: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

# Remove columns with all 0s
zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
raw_counts = raw_counts[:, ~zero_genes]
print('Counts matrix shape after removing unexpressed genes: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))


if args.filtered_barcodes is None:
    print("Will not filter barcodes as no barcode filtering file was provided")
else:
    if os.path.exists(args.filtered_barcodes):
        print("Filtered barcodes exist.")
        barcodes_filtered_df = read10x.read_barcodes(args.filtered_barcodes)
        if (any(barcodes_df['Barcode'].isin(barcodes_filtered_df['Barcode']))):
            raw_counts = raw_counts[barcodes_df['Barcode'].isin(barcodes_filtered_df['Barcode'])]
            print('\nThe original number of barcodes in the counts matrix: {}. \nThe number of barcodes in the user-provided barcode list: {}.\nThe number of barcodes after filtering for user-provided barcodes: {}'.format(barcodes_df.shape[0], barcodes_filtered_df.shape[0], raw_counts.shape[0]))
        else:
            print("There are no barcodes remaining in your dataframe after filtering on the provided --filter_barcodes file.\n\
            This is what the top your original barcodes looks like:\n {} \
            \n\nAnd this is what the filtering barcodes look like:\n {} \
            Please check that the provided filter barcode file is accurate for this data and has the same format.".format(barcodes_df['Barcode'].head, barcodes_filtered_df['Barcode'].head))
            exit()
    else:
        print("Cannot read filtered barcode file, please check the directory path and try again.\nInterpreted path for filtered barcodes file: " + args.filtered_barcodes)
        exit()


clf = doubletdetection.BoostClassifier(n_iters=args.n_iterations, use_phenograph=pheno, standard_scaling=standard_scaling, verbose = True)
doublets = clf.fit(raw_counts).predict(p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)

results = pd.Series(doublets, name="DoubletDetection_DropletType")
dataframe = pd.concat([barcodes_df, results], axis=1)


dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(1.0, "doublet")
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(0.0, "singlet")

print("Writing results to {}.".format(os.path.join(args.outdir,'DoubletDetection_doublets_singlets.tsv')))

dataframe.to_csv(os.path.join(args.outdir,'DoubletDetection_doublets_singlets.tsv'), sep = "\t", index = False)


### Figures ###
doubletdetection.plot.convergence(clf, save=os.path.join(args.outdir,'convergence_test.pdf'), show=False, p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)

f3 = doubletdetection.plot.threshold(clf, save=os.path.join(args.outdir,'threshold_test.pdf'), show=False, p_step=6)


### Make summary of singlets and doublets and write to file ###
summary = pd.DataFrame(dataframe.DoubletDetection_DropletType.value_counts())
summary.index.name = 'Classification'
summary.reset_index(inplace=True)
summary = summary.rename({'DoubletDetection_DropletType': 'Droplet N'}, axis=1)

print("Writing summary to {}.".format(os.path.join(args.outdir,'DoubletDetection_summary.tsv')))
summary.to_csv(os.path.join(args.outdir,'DoubletDetection_summary.tsv'), sep = "\t", index = False)


