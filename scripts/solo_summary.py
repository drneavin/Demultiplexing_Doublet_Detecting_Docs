#!/usr/bin/env python

##### Set up the packages #####
import argparse
import pandas as pd
import numpy
import sys
import os

sys.path.append('/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/scripts/ONEK1K/hg38/refSNVs/sceQTL-Gen-Demultiplex/mods') 
import read10x

##### Parse the variables passed to python #####
parser = argparse.ArgumentParser(
    description="wrapper for scrublet for doublet detection of transcriptomic data.")
parser.add_argument("-b", "--barcodes", required = True, help = "barcodes.tsv or barcodes.tsv.gz from cellranger")
parser.add_argument("-s", "--solo_output", required = True, help = "solo output directory")
parser.add_argument("-d", "--solo_doublet_file", required = False, default = "is_doublet.npy", help = "Which of the solo doublet classification files to use; options are is_doublet.npy (default) or preds.npy (forces doublet number if provided -e in solo command)")
args = parser.parse_args()

##### Read in the barcodes and the solo results #####
barcodes_df = read10x.read_barcodes(args.barcodes)
doublet_list = numpy.load(args.solo_output + "/" + args.solo_doublet_file)
scores_list = numpy.load(args.solo_output + "/logit_scores.npy")

##### make a final dataframe of results + barcodes #####
dataframe = barcodes_df
dataframe["solo_DropletType"] = doublet_list
dataframe["solo_DropletScore"] = scores_list

##### Replace True and False with singlet and doublet #####
dataframe.solo_DropletType = dataframe.solo_DropletType.replace(True, "doublet")
dataframe.solo_DropletType = dataframe.solo_DropletType.replace(False, "singlet")

##### Write results
dataframe.to_csv(os.path.join(args.solo_output,'solo_results.tsv'), sep = "\t", index = False)



### Make summary of singlets and doublets and write to file ###
summary = pd.DataFrame(dataframe.solo_DropletType.value_counts())
summary.index.name = 'Classification'
summary.reset_index(inplace=True)
summary = summary.rename({'solo_DropletType': 'Droplet N'}, axis=1)

summary.to_csv(os.path.join(args.solo_output,'solo_summary.tsv'), sep = "\t", index = False)
