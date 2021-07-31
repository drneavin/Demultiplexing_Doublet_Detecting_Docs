Demuxlet Tutorial
===========================
Demuxlet is a genotype demultiplexing software that requires reference genotypes to be available for each individual in the pool. 
Therefore, if you don't have reference genotypes, you may want to demultiplex with one of the softwares that do not require reference genotype data
(Freemuxlet, scSplit, Souporcell or Vireo)


Data
----
This is the data that you will need to have preparede to run Demuxlet

- Reference SNP genotypes for each individual

  + We recommend imputing your SNP genotypes and filtering to exon loctions and >= 5% minor allele frequency (MAF)

- Barcode file

- Bam file (aligned reads)


.. admonition:: Optional

    - A text file with the individual ids (separated by line) as they appear in the vcf file


Run Popscle
-----------
Poscle Pileup
^^^^^^^^^^^^^
First we will need to identify the number of reads from each allele at each SNP location::

	singularity exec image.sif popscle dsc-pileup --sam <bam> --vcf <snp_genotype_vcf> --group-list <barcodes> --out <output_name>


Popscle Demuxlet
^^^^^^^^^^^^^^^^
Once you have run ``popscle pileup``, you can demultiplex your samples::

    singularity exec image.sif
