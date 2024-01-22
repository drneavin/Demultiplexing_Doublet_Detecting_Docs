
Overview of Demultiplexing Softwares
===========================================

Demultiplexing softwares use the inherent genetic differences between donors multiplexed in a single pool to assign droplets to each donor and to identify doublets.
There are five demultiplexing softwares that have different capabilities and advantages depending on your dataset.
As you can see from this table, only :ref:`Demuxlet <Demuxlet-docs>` absolutely requires reference SNP genotypes for the donors multiplexed in your pool.
However, :ref:`Souporcell <Souporcell-docs>` and :ref:`Vireo <Vireo-docs>` are also capable of accomodating reference SNP genotypes as well.

+--------------------------------------+------------------------------------------+------------------------------------------+------------------------------------------+
| Demultiplexing Software              | .. centered:: Requires Reference         | .. centered:: Can Use Reference          | .. centered:: Estimates Ambient RNA      |
|                                      | .. centered:: SNP Genotypes              | .. centered:: SNP Genotypes              |                                          |
+======================================+==========================================+==========================================+==========================================+
|:ref:`Demuxalot <Demuxalot-docs>`     | .. centered:: |:heavy_check_mark:|       | .. centered:: |:heavy_check_mark:|       | .. centered:: |:heavy_multiplication_x:| |
+--------------------------------------+------------------------------------------+------------------------------------------+------------------------------------------+
|:ref:`Demuxlet <Demuxlet-docs>`       | .. centered:: |:heavy_check_mark:|       | .. centered:: |:heavy_check_mark:|       | .. centered:: |:heavy_multiplication_x:| |
+--------------------------------------+------------------------------------------+------------------------------------------+------------------------------------------+
|:ref:`Dropulation <Dropulation-docs>` | .. centered:: |:heavy_check_mark:|       | .. centered:: |:heavy_check_mark:|       | .. centered:: |:heavy_multiplication_x:| |
+--------------------------------------+------------------------------------------+------------------------------------------+------------------------------------------+
|:ref:`Femuxlet <Freemuxlet-docs>`     | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_multiplication_x:| |
+--------------------------------------+------------------------------------------+------------------------------------------+------------------------------------------+
|:ref:`scSplit <scSplit-docs>`         | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_multiplication_x:| |
+--------------------------------------+------------------------------------------+------------------------------------------+------------------------------------------+
|:ref:`Souporcell <Souporcell-docs>`   | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_check_mark:|       | .. centered:: |:heavy_check_mark:|       |
+--------------------------------------+------------------------------------------+------------------------------------------+------------------------------------------+
|:ref:`Vireo <Vireo-docs>`             | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_check_mark:|       | .. centered:: |:heavy_check_mark:|       |
+--------------------------------------+------------------------------------------+------------------------------------------+------------------------------------------+

We highly recommend using :ref:`Souporcell <Souporcell-docs>` if only to estimate the percentage of ambient RNA in your pool.
As far as we are aware, this is the only software that leverages SNP genotype data to estimate ambient RNA in multiplexed pools and it is helpful to identify high ambient RNA which is sometimes undetectable with basic QC metrics.
We view this as supplementary to other ambient RNA methods that use the transcriptional profile to estimate and remove ambient RNA per droplet.

If you don't know which demultiplexing software(s) to run, take a look at our :ref:`Software Selection Recommendations <SoftwareSelection-docs>` based on your dataset or use our **add widget link here**