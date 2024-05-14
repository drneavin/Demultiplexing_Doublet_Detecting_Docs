#!/usr/bin/env python



import argparse
import subprocess



parser = argparse.ArgumentParser(
    description="wrapper for execution of demuxlet and cellsnp pileup.")
parser.add_argument("-s", "--samFile", required = True, help = "Indexed sam/bam file(s), comma separated multiple samples. Mode 1a & 2a: one sam/bam file with single cell. Mode 1b & 2b: one or multiple bulk sam/bam files.")
parser.add_argument("-S", "--samFileList", required = True,  help = "A list file containing bam files, each per line, for Mode 1b & 2b.")
parser.add_argument("-O", "--outDir", required = True, help = "Output directory for VCF and sparse matrices.")
parser.add_argument("-R", "--regionsVCF", required = False, default = 'CB', help = "A vcf file listing all candidate SNPs, for fetch each variants.")
parser.add_argument("-T", "--targetsVCF", required = False, default = 'UB', help = "Similar as -R, but the next position is accessed by streaming rather than indexing/jumping (like -T in samtools/bcftools mpileup).")
parser.add_argument("-b", "--barcodeFile", required = False, help = "A plain file listing all effective cell barcode.")
parser.add_argument("-i", "--sampleList", required = False, default = None, help = "A list file containing sample IDs, each per line.")
parser.add_argument("-I", "--sampleIDs", required = False, default = None, help = "Comma separated sample ids.")
parser.add_argument("--genotype", required = False, default = False, type = bool, help = "If use, do genotyping in addition to counting. True or False. Default: False.")
parser.add_argument("--gzip", required = False, default = True, type = bool, help = "If use, the output files will be zipped into BGZF format. True or False. Default: True")
parser.add_argument("--printSkipSNPs", required = False, default = False, type = bool, help = "If use, the SNPs skipped when loading VCF will be printed. True or False, default: False.")
parser.add_argument("-p", "--nproc", required = False, default = 1, help = " Number of subprocesses [1].")
parser.add_argument("--chrom", required = False, help = "The chromosomes to use, comma separated [1 to 22].")
parser.add_argument("--cellTAG", required = False, default = 'CB', help = "Tag for cell barcodes, turn off with None [CB].")
parser.add_argument("--UMItag", required = False, default = 'Auto', help = "Tag for UMI: UR, Auto, None. For Auto mode, use UR if barcodes is inputted, otherwise use None. None mode means no UMI but read counts [Auto].")
parser.add_argument("--minCOUNT", required = False, default = 20, type = int, help = "Minimum aggragated count [20].")
parser.add_argument("--minMAF", required = False, default = 0, type = float, help = "Minimum minor allele frequency [0.00].")
parser.add_argument("--doubletGL", required = False, help = "If use, keep doublet GT likelihood, i.e., GT=0.5 and GT=1.5.")
parser.add_argument("--inclFLAG", required = False, help = "Required flags: skip reads with all mask bits unset []")
parser.add_argument("--exclFLAG", required = False, default = 'UNMAP,SECONDARY,QCFAIL', help = "Filter flags: skip reads with any mask bits set [UNMAP,SECONDARY,QCFAIL (when use UMI) or UNMAP,SECONDARY,QCFAIL,DUP (otherwise)].")
parser.add_argument("--minLEN", required = False, default = 30, help = "Minimum mapped length for read filtering [30].")
parser.add_argument("--minMAPQ", required = False, default = 20, help = "Minimum MAPQ for read filtering [20].")
parser.add_argument("--countORPHAN", required = False, default = False, help = "If use, do not skip anomalous read pairs. True or False. Default: False.")
args = parser.parse_args()

print("checking genome nomenclature for vcf and bam files.")

return_code = subprocess.run("compare_vcf_bam_genome.sh " + args.bam + " " + args.vcf, shell = True).returncode

if return_code != 0:
    exit(0)




if args.genotype:
    genotype = ' --genotype '
else:
    genotype = ''
    
if args.gzip:
    gzip = ' --gzip '
else:
    gzip = ''

if args.printSkipSNPs:
    printSkipSNPs = ' --printSkipSNPs '
else:
    printSkipSNPs = ''

if args.countORPHAN:
    countORPHAN = ' --countORPHAN '
else:
    countORPHAN = ''

print("Running cellsnp-lite pileup.")


subprocess.run("cellsnp_pileup.py -s " + args.bam +
                " -S " + args.samFileList +
                " -O " + args.outDir +
                " -R " + args.regionsVCF +
                " -T " + args.targetsVCF +
                " -b " + args.barcodeFile +
                " -i " + args.sampleList +
                " -I " + args.sampleIDs +
                genotype +
                gzip +
                printSkipSNPs +
                " -p " + args.nproc +
                " --chrom " + args.chrom +
                " --cellTAG " + args.cellTAG +
                " --UMItag " + args.UMItag +
                " --minCOUNT " + args.minCOUNT +
                " --minMAF " + args.minMAF +
                " --doubletGL " + args.doubletGL +
                " --inclFLAG " + args.inclFLAG +
                " --exclFLAG " + args.exclFLAG +
                " --minLEN " + args.minLEN +
                " --minMAPQ " + args.minMAPQ +
                countORPHAN, shell = True)