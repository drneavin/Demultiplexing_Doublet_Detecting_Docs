#!/usr/bin/env python

import numpy as np
import demuxalot
from demuxalot import Demultiplexer, BarcodeHandler, ProbabilisticGenotypes, count_snps
import pandas as pd
import argparse

print("loaded packages")

parser = argparse.ArgumentParser(
    description="wrapper for demuxalot for demultiplexing of multiplexed single-cell data.")
parser.add_argument("-b", "--barcodes", required = True, help = "Input barcodes file")
parser.add_argument("-a", "--bamfile", required = True, help = "Input bam file")
parser.add_argument("-n", "--indiv_file", required = True, default = None, help = "File containing individuals in pool separated by line. No header.")
parser.add_argument("-v", "--vcf", required = True, default = None, help = "The vcf file that has the SNP genotypes for each donor in the pool.")
parser.add_argument("-o", "--outdir", required = True, default = None, help = "The output directory.")
parser.add_argument("-r", "--refine", required = True, default = True, help = "Whether to run genotype refinement.")
parser.add_argument("-p", "--nproc", required = False, type = int, default = 1, help = "(max.) number of parallel jobs; -1: number of CPU cores/threads; default: 1.")
args = parser.parse_args()

print("read arguments")



### Load in file with individual IDs as list
with open(args.indiv_file) as f:
    names = f.read().splitlines()

names_ordered = np.sort(names)

# Load genotypes
print("loading genotypes")
genotypes = ProbabilisticGenotypes(genotype_names=names_ordered)
genotypes.add_vcf(args.vcf)


# Load barcodes
print("loading barcodes")
barcode_handler = BarcodeHandler.from_file(args.barcodes)

snps = count_snps(
    bamfile_location = args.bamfile,
    chromosome2positions=genotypes.get_chromosome2positions(),
    barcode_handler=barcode_handler, 
    joblib_n_jobs=args.nproc,
)

print("estimating posterior probs and likelihoods")
# returns two dataframes with likelihoods and posterior probabilities 
likelihoods, posterior_probabilities = Demultiplexer.predict_posteriors(
    snps,
    genotypes=genotypes,
    barcode_handler=barcode_handler,
)

print("writing unrefined results")
# Write output
likelihoods.to_csv(args.outdir + "/likelihoods.tsv.gz", sep = "\t", index = True)
posterior_probabilities.to_csv(args.outdir + "/posterior_probabilities.tsv.gz", sep = "\t", index = True)
posterior_probabilities.idxmax(axis = 1).to_csv(args.outdir + "/assignments.tsv.gz", sep = "\t", index = True)

if args.refine:

    print("inferring refined genotypes")
    # Infer refined genotypes 
    refined_genotypes, _posterior_probabilities = Demultiplexer.learn_genotypes(snps, genotypes=genotypes, n_iterations=5, barcode_handler=barcode_handler)

    print("estimating likelihoods and probabilities")
    # Use learnt genotypes for demultiplexing
    likelihoods, posterior_probabilities = Demultiplexer.predict_posteriors(
        snps,
        genotypes=refined_genotypes,
        barcode_handler=barcode_handler,
    )

    print("writing second set of output")
    # Write output
    likelihoods.to_csv(args.outdir + "/likelihoods_refined.tsv.gz", sep = "\t", index = True)
    posterior_probabilities.to_csv(args.outdir + "/posterior_probabilities_refined.tsv.gz", sep = "\t", index = True)
    posterior_probabilities.idxmax(axis = 1).to_csv(args.outdir + "/assignments_refined.tsv.gz", sep = "\t", index = True)

print("Done!")
