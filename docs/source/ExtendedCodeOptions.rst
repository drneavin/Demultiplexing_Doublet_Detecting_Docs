Extended Code Options
=====================
The extended options for each of the softwares outlined in this workflow are provided here as a courtesy.
You can, of course, find extended details about each of the softwares at their respective websites


.. _cellSNP:

cellSNP
~~~~~~~~

.. code-block:: bash

    cellSNP -s $BAM -b $BARCODES -o $OUTDIR/vireo_cellsnp.vcf -R $VCF -p 20 --minMAF 0.1 --minCOUNT 20
