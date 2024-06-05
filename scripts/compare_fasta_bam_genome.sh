#!/bin/bash


BAM=$1
FASTA=$2


# BAM=/directflow/SCCGGroupShare/projects/data/experimental_data/RAW/BreastQTL/scRNA/Mastectomy_pool4/outs/possorted_genome_bam.bam ## yes chr
# FASTA=/directflow/SCCGGroupShare/projects/DrewNeavin/References/ENSEMBLfasta/GRCh38/genome.fa ## no chr
# /directflow/SCCGGroupShare/projects/DrewNeavin/References/UCSCrefs/hg38/hg38.fa


if $(samtools view $BAM | head -n 1 | awk '{print $3}' | grep -q 'chr')
then
    echo "BAM file is using chr notation"
    bam_genome="UCSC"
else
    echo "BAM file is using non-chr notation"
    bam_genome="ENSEMBL/NCBI"
fi



if $(head $FASTA -n 1  | awk '{print $1}' | sed 's/>//g' | grep -q "chr")
then
    echo "FASTA file is using chr notation"
    fasta_genome="UCSC"
else
    echo "FASTA file is using non-chr notation"
    fasta_genome="ENSEMBL/NCBI"
fi


if [ $bam_genome == $fasta_genome ]
then
    echo "BAM and FASTA files are using the same genome notation"
    exit 0
else
    echo "ERROR: BAM and FASTA files are using different genome notations."
    echo "Please update so that the same genome notation is used for both files."
    echo "i.e. chr1, chr2, chr3... vs 1, 2, 3."
    echo "Exiting."
    exit 1
fi