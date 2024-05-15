#!/usr/bin/env python



import argparse
import subprocess
import os

__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))


parser = argparse.ArgumentParser(
    description="wrapper for execution of dropulation demultiplexing.")
parser.add_argument("--CELL_BC_FILE", required = False, default = None, help = "Override NUM_CORE_BARCODES and process reads that have the cell barcodes in this file instead. The file has 1 column with no header. Required. Cannot be used in conjunction with argument(s) NUM_BARCODES.")
parser.add_argument("-I", "--INPUT_BAM", required = True, default = None, help = "Input bam file.")
parser.add_argument("--NUM_BARCODES", required = False,  default = None, help = "Number of cells that you want to extract from the BAM. The program will pick the top <NUM_BARCODES> barcodes by read count. Required. Cannot be used in conjunction with argument(s) CELL_BC_FILE")
parser.add_argument("-O", "--OUTPUT", required = True, help = "Output file of sample assignments. This supports zipped formats like gz and bz2. Required.")
parser.add_argument("--VCF", required = True, help = "The input VCF file to analyze. Required.")
parser.add_argument("--ADD_MISSING_VALUES", required = False, default = 'true', help = "<Boolean>on a per-snp basis, generate a mxiture of all the data across known samples to fill in likelihoods for the missing samples.  Default value: true. Possible values: {true, false}.")
parser.add_argument("--ALLELE_FREQUENCY_ESTIMATE_FILE", required = False, help = "A file that contains an estimate of the allele frequency expected for each SNP across donors. The best estimate of this will come from the fraction of reference and alternate alleleUMIs that are observed at each snp site.  This report can be generated via GatherDigitalAlleleCounts.  This is a fractional estimate between 0 and 1.  File is tab seperated, with at least 3 columns:chromosome, position, maf_umi.  When supplied and CELL_CONTAMINATION_ESTIMATE_FILE is provided, this modifies the likelihood error rates to take into account how often the allele observed can be drawn from ambient RNA.  Default value: null.")
parser.add_argument("--ANSWER_KEY_FILE", required = False, help = "This file provides an answer key for the input BAM and VCF. This does not change any of the likelihood analysis results, but appends an additional column to the output that contains the known sample assignment for a cell.This is useful when assessing how well the method works, but is completely optional. The format of the file is 2 tab seperated columns with a header [cell sample]. The first column contains the cell barcodes from the BAM, the second the sample assignments from the VCF. Default value: null..")
parser.add_argument("--arguments_file", required = False, help = "read one or more arguments files and add them to the command line  This argument may be specified 0 or more times. Default value: null. .")
parser.add_argument("--BAM_OUTPUT", required = False, help = "Output a version of the BAM containing only the reads that have coverage in the input BAM and VCF.  Default value: null.")
parser.add_argument("--CELL_BARCODE_TAG", required = False, default = 'XC', help = "The cell barcode tag. If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode. Default value: XC.")
parser.add_argument("--CELL_CONTAMINATION_ESTIMATE_FILE", required = False,  help = "A file that contains an estimate of how much ambient RNA is in each cell. This is a fractional estimate between 0 and 1. File is tab seperated, with 2 columns:cell_barcode and frac_contamination. When supplied along side the ALLELE_FREQUENCY_ESTIMATE_FILE, this modifies the likelihood error rates to take into account how oftenthe allele observed can be drawn from ambient RNA. We use cellbender remove background [https://github.com/broadinstitute/CellBender] to estimate the number of transcripts before and after ambient cleanup to define the fraction of transcripts that come from ambient RNA. Default value: null.")
parser.add_argument("--COMPRESSION_LEVEL", required = False, default = 5, help = "Compression level for all compressed files created (e.g. BAM and VCF).  Default value: 5.")
parser.add_argument("--CREATE_INDEX", required = False, default = 'false', help = "Whether to create an index when writing VCF or coordinate sorted BAM output.  Default value: false. Possible values: {true, false} .")
parser.add_argument("--CREATE_MD5_FILE", required = False, default = 'false', help = "Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. Possible values: {true, false.")
parser.add_argument("--DNA_MODE", required = False, default = 'false', help = "EXPERIMENTAL!!! Run the program in DNA Mode. In this mode, reads should have a cell barcode, but will be missing gene annotations and UMIs. All reads will be accepted as passing, and each read (or read pair) will be treated as a single UMI If the data is PCR Duplicate marked, duplicate reads will be filtered. Default value: false. Possible values: {true, false}.")
parser.add_argument("--EDIT_DISTANCE", required = False, default = 1, help = "The edit distance that molecular barcodes should be combined at within a gene/SNP. Default value: 1.")
parser.add_argument("--FIXED_ERROR_RATE", required = False, help = "Instead of using base qualities to determine error rate, use a fixed error rate instead. This is rounded to the nearest phread score internally.  Default value: null..")
parser.add_argument("--FRACTION_SAMPLES_PASSING", required = False, default = 0.5, help = "At least <FRACTION_SAMPLES_PASSING> samples must have genotype scores >= GQ_THRESHOLD for the variant in the VCF to be included in the analysis.  Default value: 0.5.")
parser.add_argument("--FUNCTION_TAG", required = False, default = 'XF', help = "The functional annotation for the read. If set, extracts the functional annotation(s) [CODING/UTR/etc] at each SNP position and outputs in the verbose output. Default value: XF.")
parser.add_argument("--GA4GH_CLIENT_SECRETS", required = False, default = 0, help = "Google Genomics API client_secrets.json file path.  Default value: client_secrets.json.")
parser.add_argument("--GENE_FUNCTION_TAG", required = False, default = 'gf', help = "Gene Function tag.  For a given gene name <GENE_NAME_TAG>, this is the function of the gene at this read's position: UTR/CODING/INTRONIC/...  Default value: gf.")
parser.add_argument("--GENE_NAME_TAG", required = False, default = 'gn', help = "Gene Name tag.  Takes on the gene name this read overlaps (if any)  Default value: gn.")
parser.add_argument("--GENE_STRAND_TAG", required = False, default = 'gs', help = "Gene Strand tag.  For a given gene name <GENE_NAME_TAG>, this is the strand of the gene. Default value: gs.")
parser.add_argument("--GQ_THRESHOLD", required = False, default = 30, help = "The minimum genotype quality for a variant.  Set this value to 0 to not filter by GQ scores if they are present, or to -1 to completely ignore GQ values if they are not set in the genotype info field. If the GQ field is not set in the VCF header, this will be set to -1 by default. Default value: 30.")
parser.add_argument("--IGNORED_CHROMOSOMES", required = False,  help = "A list of chromosomes to omit from the analysis. The default is to omit the sex chromosomes. This argument may be specified 0 or more times. Default value: [X, Y, MT].")
parser.add_argument("--LOCUS_FUNCTION_LIST", required = False, help = "A list of functional annotations that reads need to be completely contained by to be considered for analysis. This argument may be specified 0 or more times. Default value: [CODING, UTR]. Possible values: {INTERGENIC, INTRONIC, UTR, CODING, RIBOSOMAL}.")
parser.add_argument("--MAX_ERROR_RATE", required = False, help = "Caps the base error rate at a maximum probability so no SNP can be weighed more than this value. For example, if this value was 0.01, then a base quality 30 value (normally an erro rate of 0.001) would become 0.01. With the same threshold, a base with an error rate of 0.1 would be unaffected. Default value: null.")
parser.add_argument("--MAX_RECORDS_IN_RAM", required = False, default = 500000, help = "When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.  Default value: 500000.")
parser.add_argument("--MOLECULAR_BARCODE_TAG", required = False, default = 'SM', help = "The molecular barcode tag.  Default value: XM.")
parser.add_argument("--QUIET", required = False, default = 'false', help = "Whether to suppress job-summary info on System.err.  Default value: false. Possible values: {true, false}.")
parser.add_argument("--READ_MQ", required = False, default = 10, help = "The map quality of the read to be included.  Default value: 10.")
parser.add_argument("--REFERENCE_SEQUENCE", required = False,  help = "Reference sequence file.  Default value: null.")
parser.add_argument("--RETAIN_MONOMORPIC_SNPS", required = False, default = 'false', help = "Should monomorphic SNPs across the population be retained?  Default value: false. Possible values: {true, false}.")
parser.add_argument("--SAMPLE_FILE", required = False, help = "A file with a list of samples in the VCF to match up to cells in the BAM. This subsets the VCF into a smaller data set containing only the samples listed. The file has 1 column with no header.  Default value: null.")
parser.add_argument("--SNP_LOG_RATE", required = False, default = 1000, help = "Logging interval for SNP sequence pileup processing  Default value: 1000.")
parser.add_argument("--STRAND_STRATEGY", required = False, default = 'SENSE', help = "The strand strategy decides which reads will be used by analysis. The SENSE strategy requires the read and annotation to have the same strand. The ANTISENSE strategy requires the read and annotation to be on opposite strands. The BOTH strategy is permissive, and allows the read to be on either strand. Default value: SENSE. Possible values: {SENSE, ANTISENSE, BOTH}.")
parser.add_argument("--TMP_DIR", required = False, help = "One or more directories with space available to be used by this program for temporary storage of working files. This argument may be specified 0 or more times. Default value: null.")
parser.add_argument("--TRANSCRIBED_SNP_FAIL_FAST_THRESHOLD", required = False, default = -1, help = "(Optional) If set to a value, this program will quit early with an exception if this number of UMIs are observed without encountering a transcribed SNP - that is, a variant from the VCF file that is will be informative during donor assignment. A value of -1 disables this test. As some fraction of UMIs may notcontain transcribed SNPs, the value for this should probably be set to a fairly large number, like 10% of the total UMIs in your experiment. Mostly useful for debugging / or testing new data sets. Default value: -1.")
parser.add_argument("--USE_JDK_DEFLATER", required = False, default = 'false', help = "Use the JDK Deflater instead of the Intel Deflater for writing compressed output. Default value: false. Possible values: {true, false}.")
parser.add_argument("--USE_JDK_INFLATER", required = False, default = 'false', help = "Use the JDK Inflater instead of the Intel Inflater for reading compressed input. Default value: false. Possible values: {true, false}.")
parser.add_argument("--VALIDATION_STRINGENCY", required = False, default = 'STRICT', help = "Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT}.")
parser.add_argument("--VCF_OUTPUT", required = False, help = "Output a version of the VCF containing only the variants that have coverage in the input BAM. All samples are retained. This file can be used as the VCF input file on subsequent runs of the same data set to greatly speed up the run, if the same BAM is being analyzed with different conditions. If you plan on calling doublets with DetectDoublets, you'll need this file. Default value: null.")
parser.add_argument("--VERBOSE_BEST_DONOR_OUTPUT", required = False, help = "Verbose output of every cell/SNP result for the BEST donor. This supports zipped formats like gz and bz2. Default value: null.")
parser.add_argument("--VERBOSE_OUTPUT", required = False, help = "Verbose output of every cell/SNP result. This supports zipped formats like gz and bz2. Default value: null.")
parser.add_argument("--VERBOSITY", required = False, default = 'INFO', help = "display the version number for this tool  Default value: true. Allowed values are Possible values: {ERROR, WARNING, INFO, DEBUG}.")
args = parser.parse_args()


print()


print("checking genome nomenclature for vcf and bam files.")

return_code = subprocess.run(os.path.join(__location__, "compare_vcf_bam_genome.sh ") + args.INPUT_BAM + " " + args.VCF, shell = True).returncode

if return_code != 0:
    exit(0)

if args.CELL_BC_FILE is None:
    CELL_BC_FILE = ''
    if args.NUM_BARCODES is None:
        print("Error: Either CELL_BC_FILE or NUM_BARCODES must be specified.")
        exit(0)
else:
    CELL_BC_FILE = (' --CELL_BC_FILE ' + str(args.CELL_BC_FILE))

if args.NUM_BARCODES is None:
    NUM_BARCODES = ''
else:
    NUM_BARCODES = (' --NUM_BARCODES ' + str(args.NUM_BARCODES))

## Check for optional arguments
if args.ALLELE_FREQUENCY_ESTIMATE_FILE is None:
    ALLELE_FREQUENCY_ESTIMATE_FILE = ''
else:
    ALLELE_FREQUENCY_ESTIMATE_FILE = (' --ALLELE_FREQUENCY_ESTIMATE_FILE ' + str(args.ALLELE_FREQUENCY_ESTIMATE_FILE))

if args.ANSWER_KEY_FILE is None:
    ANSWER_KEY_FILE = ''
else:
    ANSWER_KEY_FILE = (' --ANSWER_KEY_FILE ' + str(args.ANSWER_KEY_FILE))

if args.arguments_file is None:
    arguments_file = ''
else:   
    arguments_file = (' --arguments_file ' + str(args.arguments_file))

if args.BAM_OUTPUT is None:
    BAM_OUTPUT = ''
else:
    BAM_OUTPUT = (' --BAM_OUTPUT ' + str(args.BAM_OUTPUT))

if args.CELL_CONTAMINATION_ESTIMATE_FILE is None:
    CELL_CONTAMINATION_ESTIMATE_FILE = ''
else:
    CELL_CONTAMINATION_ESTIMATE_FILE = (' --CELL_CONTAMINATION_ESTIMATE_FILE ' + str(args.CELL_CONTAMINATION_ESTIMATE_FILE))

if args.FIXED_ERROR_RATE is None:
    FIXED_ERROR_RATE = ''  
else:
    FIXED_ERROR_RATE = (' --COMPRESSION_LEVEL ' + str(args.FIXED_ERROR_RATE))

if args.MAX_ERROR_RATE is None:
    MAX_ERROR_RATE = ''
else:
    MAX_ERROR_RATE = (' --MAX_ERROR_RATE ' + str(args.MAX_ERROR_RATE))

if args.REFERENCE_SEQUENCE is None:
    REFERENCE_SEQUENCE = ''
else:
    REFERENCE_SEQUENCE = (' --REFERENCE_SEQUENCE ' + str(args.REFERENCE_SEQUENCE))

if args.SAMPLE_FILE is None:
    SAMPLE_FILE = ''
else:
    SAMPLE_FILE = (' --SAMPLE_FILE ' + str(args.SAMPLE_FILE))

if args.TMP_DIR is None:
    TMP_DIR = ''
else:
    TMP_DIR = (' --TMP_DIR ' + str(args.TMP_DIR))
    
if args.VCF_OUTPUT is None:
    VCF_OUTPUT = ''
else:
    VCF_OUTPUT = (' --VCF_OUTPUT ' + str(args.VCF_OUTPUT))

if args.VERBOSE_BEST_DONOR_OUTPUT is None:
    VERBOSE_BEST_DONOR_OUTPUT = ''
else:
    VERBOSE_BEST_DONOR_OUTPUT = (' --VERBOSE_BEST_DONOR_OUTPUT ' + str(args.VERBOSE_BEST_DONOR_OUTPUT))

if args.VERBOSE_OUTPUT is None:
    VERBOSE_OUTPUT = ''
else:
    VERBOSE_OUTPUT = (' --VERBOSE_OUTPUT ' + str(args.VERBOSE_OUTPUT))

if args.LOCUS_FUNCTION_LIST is None:
    LOCUS_FUNCTION_LIST = ''
else:
    LOCUS_FUNCTION_LIST = (' --LOCUS_FUNCTION_LIST ' + str(args.LOCUS_FUNCTION_LIST))

if args.IGNORED_CHROMOSOMES is None:
    IGNORED_CHROMOSOMES = ''
else:
    IGNORED_CHROMOSOMES = (' --IGNORED_CHROMOSOMES ' + str(args.IGNORED_CHROMOSOMES))
    
print(CELL_BC_FILE)

subprocess.run("AssignCellsToSamples " + CELL_BC_FILE + 
                " -I " + str(args.INPUT_BAM) + 
                " -O " + str(args.OUTPUT) + 
                NUM_BARCODES +
                " --VCF " + str(args.VCF) + 
                " --ADD_MISSING_VALUES " + str(args.ADD_MISSING_VALUES) +
                ALLELE_FREQUENCY_ESTIMATE_FILE +
                ANSWER_KEY_FILE +
                arguments_file +
                BAM_OUTPUT +
                " --CELL_BARCODE_TAG " + str(args.CELL_BARCODE_TAG) +
                CELL_CONTAMINATION_ESTIMATE_FILE +
                " --COMPRESSION_LEVEL " + str(args.COMPRESSION_LEVEL) +
                " --CREATE_INDEX " + str(args.CREATE_INDEX) +
                " --CREATE_MD5_FILE " + str(args.CREATE_MD5_FILE) +
                " --DNA_MODE " + str(args.DNA_MODE) +
                " --EDIT_DISTANCE " + str(args.EDIT_DISTANCE) +
                FIXED_ERROR_RATE +
                " --FRACTION_SAMPLES_PASSING " + str(args.FRACTION_SAMPLES_PASSING) +
                " --FUNCTION_TAG " + str(args.FUNCTION_TAG) +
                " --GA4GH_CLIENT_SECRETS " + str(args.GA4GH_CLIENT_SECRETS) +
                " --GENE_FUNCTION_TAG " + str(args.GENE_FUNCTION_TAG) +
                " --GENE_NAME_TAG " + str(args.GENE_NAME_TAG) +
                " --GENE_STRAND_TAG " + str(args.GENE_STRAND_TAG) +
                " --GQ_THRESHOLD " + str(args.GQ_THRESHOLD) +
                IGNORED_CHROMOSOMES +
                LOCUS_FUNCTION_LIST +
                MAX_ERROR_RATE +
                " --MAX_RECORDS_IN_RAM " + str(args.MAX_RECORDS_IN_RAM) +
                " --MOLECULAR_BARCODE_TAG " + str(args.MOLECULAR_BARCODE_TAG) +
                " --QUIET " + str(args.QUIET) +
                " --READ_MQ " + str(args.READ_MQ) +
                REFERENCE_SEQUENCE +
                " --RETAIN_MONOMORPIC_SNPS " + str(args.RETAIN_MONOMORPIC_SNPS) +
                SAMPLE_FILE +
                " --SNP_LOG_RATE " + str(args.SNP_LOG_RATE) +
                " --STRAND_STRATEGY " + str(args.STRAND_STRATEGY) +
                TMP_DIR +
                " --TRANSCRIBED_SNP_FAIL_FAST_THRESHOLD " + str(args.TRANSCRIBED_SNP_FAIL_FAST_THRESHOLD) +
                " --USE_JDK_DEFLATER " + str(args.USE_JDK_DEFLATER) +
                " --USE_JDK_INFLATER " + str(args.USE_JDK_INFLATER) +
                " --VALIDATION_STRINGENCY " + str(args.VALIDATION_STRINGENCY) +
                VCF_OUTPUT +
                VERBOSE_BEST_DONOR_OUTPUT +
                VERBOSE_OUTPUT +
                " --VERBOSITY " + str(args.VERBOSITY), 
                shell = True)
