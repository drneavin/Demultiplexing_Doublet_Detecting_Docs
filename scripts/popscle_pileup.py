#!/usr/bin/env python


import argparse
import subprocess
import os

__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

parser = argparse.ArgumentParser(
    description="wrapper for execution of demuxlet and freemuxlet popscle pileup.")
parser.add_argument("--sam", required = True, help = "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed.")
parser.add_argument("--vcf", required = True,  help = "Input VCF/BCF file, containing the AC and AN field.")
parser.add_argument("--out", required = True, help = "Output file prefix.")
parser.add_argument("--tag-group", required = False, default = 'CB', help = "File containing a filtered list of droplet barcodes. This may be used if you want to use a filtered list of barcodes for doublet detection (ie need to remove droplets that are empty or high in ambient RNA).")
parser.add_argument("--tag-UMI", required = False, default = 'UB', help = " Tag representing UMIs. For 10x genomiucs, use UB.")
parser.add_argument("--exclude-flag", required = False, type = int, default = 1796, help = "SAM/BAM flag to exclude.")
parser.add_argument("--sm", required = False, default = None, help = "List of sample IDs to compare to (default: use all).")
parser.add_argument("--sm-list", required = False, default = None, help = "File containing the list of sample IDs to compare.")
parser.add_argument("--sam-verbose", required = False, default = 1000000, type = int, help = "Verbose message frequency for SAM/BAM/CRAM.")
parser.add_argument("--vcf-verbose", required = False, default = 10000, type = int, help = "Verbose message frequency for VCF/BCF.")
parser.add_argument("--skip-umi", required = False, default = False, type = bool, help = "Do not generate [prefix].umi.gz file, which stores the regions covered by each barcode/UMI pair. True or False. Default: False")
parser.add_argument("--cap-BQ", required = False, default = 40, help = "Maximum base quality (higher BQ will be capped).")
parser.add_argument("--min-BQ", required = False, default = 13, help = "Minimum base quality to consider (lower BQ will be skipped).")
parser.add_argument("--min-MQ", required = False, default = 20, help = "Minimum mapping quality to consider (lower MQ will be ignored).")
parser.add_argument("--min-TD", required = False, default = 0, help = "Minimum distance to the tail (lower will be ignored).")
parser.add_argument("--excl-flag", required = False, default = 3844, help = "SAM/BAM FLAGs to be excluded.")
parser.add_argument("--group-list", required = False, default = None, help = "List of tag readgroup/cell barcode to consider in this run. All other barcodes will be ignored. This is useful for parallelized run.")
parser.add_argument("--min-total", required = False, default = 0, help = "Minimum number of total reads for a droplet/cell to be considered.")
parser.add_argument("--min-uniq", required = False, default = 0, help = "Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered.")
parser.add_argument("--min-snp", required = False, default = 0, help = "Minimum number of SNPs with coverage for a droplet/cell to be considered.")
args = parser.parse_args()

print(args)

print("checking genome nomenclature for vcf and bam files.")

return_code = subprocess.run(os.path.join(__location__, "compare_vcf_bam_genome.sh ") + args.sam + " " + args.vcf, shell = True).returncode

if return_code != 0:
    exit(0)


##### deal with the optional arguments
if args.sm is None:
    sm = ''
else:
    sm = ' --sm ' + args.sm
           
if args.sm_list is None:
    sm_list = ''
else:
    sm_list = ' --sm-list ' + args.sm_list

if args.group_list is None:
    group_list = ''
else:
    group_list = ' --group-list ' + args.group_list

if args.skip_umi is None:
    skip_umi = ''
else:
    skip_umi = ' --skip-umi '


print("Running popscle pileup.")

subprocess.run("popscle dsc-pileup --sam " + args.sam + 
               " --vcf " + args.vcf + 
               " --out " + args.out + 
               " --tag-group " + args.tag_group + 
               " --tag-UMI " + args.tag_UMI + 
               " --exclude-flag " + str(args.exclude_flag) + 
               sm + 
               sm_list +
               " --sam-verbose " + str(args.sam_verbose) + 
               " --vcf-verbose " + str(args.vcf_verbose) + 
               skip_umi + 
               " --cap-BQ " + str(args.cap_BQ) + 
               " --min-BQ " + str(args.min_BQ) + 
               " --min-MQ " + str(args.min_MQ) + 
               " --min-TD " + str(args.min_TD) + 
               " --excl-flag " + str(args.excl_flag) + 
               group_list + 
               " --min-total " + str(args.min_total) + 
               " --min-uniq " + str(args.min_uniq) + 
               " --min-snp " + str(args.min_snp), shell = True)

