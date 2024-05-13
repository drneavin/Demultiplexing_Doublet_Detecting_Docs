#!/bin/bash


BAM=$1
VCF=$2


# BAM=/directflow/SCCGGroupShare/projects/data/experimental_data/RAW/BreastQTL/scRNA/Mastectomy_pool4/outs/possorted_genome_bam.bam ## yes chr
# VCF=/directflow/SCCGGroupShare/projects/SNP_Genotype_Processing/Imputation/12_Powell_250324/data/Imputation_input.vcf ## no chr


if $(samtools view $BAM | head -n 1 | awk '{print $3}' | grep -q 'chr')
then
    echo "BAM file is using chr notation"
    bam_genome="UCSC"
else
    echo "BAM file is using non-chr notation"
    bam_genome="ENSEMBL/NCBI"
fi



if $(grep -v "#" $VCF | head -n 1 | grep -q "^chr")
then
    echo "VCF file is using chr notation"
    vcf_genome="UCSC"
else
    echo "VCF file is using non-chr notation"
    vcf_genome="ENSEMBL/NCBI"
fi


if [ $bam_genome == $vcf_genome ]
then
    echo "BAM and VCF files are using the same genome notation"
    exit 0
else
    echo "ERROR: BAM and VCF files are using different genome notations."
    echo "Please update so that the same genome notation is used for both files."
    echo "i.e. chr1, chr2, chr3... vs 1, 2, 3."
    echo "Exiting."
    exit 1
fi