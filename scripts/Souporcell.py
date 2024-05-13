#!/usr/bin/env python


import argparse
import subprocess

parser = argparse.ArgumentParser(
    description="wrapper for execution of demuxlet and souporcell.")
parser.add_argument("-i", "--bam", required = True, help = "cellranger bam")
parser.add_argument("-b", "--barcodes", required = True, help = "barcodes.tsv from cellranger")
parser.add_argument("-f", "--fasta", required = True, help = "reference fasta file")
parser.add_argument("-t", "--threads", required = True, type = int, help = "max threads to use")
parser.add_argument("-o", "--out_dir", required = True, help = "name of directory to place souporcell files")
parser.add_argument("-k", "--clusters", required = True, help = "number cluster, tbd add easy way to run on a range of k")
parser.add_argument("-p", "--ploidy", required = False, default = "2", help = "ploidy, must be 1 or 2, default = 2")
parser.add_argument("--min_alt", required = False, default = "10", help = "min alt to use locus, default = 10.")
parser.add_argument("--min_ref", required = False, default = "10", help = "min ref to use locus, default = 10.")
parser.add_argument("--max_loci", required = False, default = "2048", help = "max loci per cell, affects speed, default = 2048.")
parser.add_argument("--restarts", required = False, default = 100, type = int, 
    help = "number of restarts in clustering, when there are > 12 clusters we recommend increasing this to avoid local minima")
parser.add_argument("--common_variants", required = False, default = None, 
    help = "common variant loci or known variant loci vcf, must be vs same reference fasta")
parser.add_argument("--known_genotypes", required = False, default = None, 
    help = "known variants per clone in population vcf mode, must be .vcf right now we dont accept gzip or bcf sorry")
parser.add_argument("--known_genotypes_sample_names", required = False, nargs = '+', default = None, 
    help = "which samples in population vcf from known genotypes option represent the donors in your sample")
parser.add_argument("--skip_remap", required = False, default = False, type = bool, 
    help = "don't remap with minimap2 (not recommended unless in conjunction with --common_variants")
parser.add_argument("--no_umi", required = False, default = "False", help = "set to True if your bam has no UMI tag, will ignore/override --umi_tag")
parser.add_argument("--umi_tag", required = False, default = "UB", help = "set if your umi tag is not UB")
parser.add_argument("--cell_tag", required = False, default = "CB", help = "DOES NOT WORK, vartrix doesnt support this! set if your cell barcode tag is not CB")
parser.add_argument("--ignore", required = False, default = False, type = bool, help = "set to True to ignore data error assertions")
parser.add_argument("--aligner", required = False, default = "minimap2", help = "optionally change to HISAT2 if you have it installed, not included in singularity build")
args = parser.parse_args()



print("checking genome nomenclature for vcf and bam files.")

return_code = subprocess.run("/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Demultiplexing_Doublet_Detecting_Docs/scripts/compare_vcf_bam_genome.sh " + args.bam + " " + args.vcf, shell = True).returncode
return_code2 = subprocess.run("/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Demultiplexing_Doublet_Detecting_Docs/scripts/compare_fasta_bam_genome.sh " + args.bam + " " + args.fasta, shell = True).returncode

if return_code != 0:
    exit(0)


subprocess.run("souporcell_pipeline.py --bam " + args.bam + 
                " --barcodes " + args.barcodes + 
                " --fasta " + args.fasta + 
                " --threads " + str(args.threads) + 
                " --out_dir " + args.out_dir + 
                " --clusters " + args.clusters + 
                " --ploidy " + args.ploidy + 
                " --min_alt " + args.min_alt + 
                " --min_ref " + args.min_ref + 
                " --max_loci " + args.max_loci + 
                " --restarts " + str(args.restarts) + 
                " --common_variants " + args.common_variants + 
                " --known_genotypes " + args.known_genotypes + 
                " --known_genotypes_sample_names " + args.known_genotypes_sample_names +
                 " --skip_remap " + str(args.skip_remap) + 
                 " --no_umi " + str(args.no_umi) + 
                 " --umi_tag " + args.umi_tag + 
                 " --cell_tag " + args.cell_tag + 
                 " --ignore " + str(args.ignore) + 
                 " --aligner " + args.aligner, shell = True)

print("souporcell pipeline complete.")