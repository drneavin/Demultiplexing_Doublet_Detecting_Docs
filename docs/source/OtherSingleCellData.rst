Considerations for Other Single Cell Data Types
===============================================

This workflow was designed for demultiplexing and detecting doublets in scRNA-seq data.
However, additional data types are becoming more frequently used - `i.e.` snRNA-seq, snATAC-seq and dual snRNA-seq + scATAC-seq.
Based on our experiences with this data we have some recommendations and considerations to take into account when applying demultiplexing and doublet detecting softwares to these data types.


.. _snrna:

snRNA-seq
---------
Demultiplexing Softwares
^^^^^^^^^^^^^^^^^^^^^^^^
We have not tested any demultiplexing softwares on snRNA-seq data in our hands but we anticipate that it should behave similarly to scRNA-seq.
The only difference we would suggest is to filter SNPs overlapping **genes** instead of just overlapping **exons**.
If you are running in to any issues or would like a discussion about use of demultiplexing softwares on snRNA-seq data, please feel free to reach out.


Doublet Detecting Softwares
^^^^^^^^^^^^^^^^^^^^^^^^^^^
We have not tested doublet detecting softwares on snRNA-seq data but the softwares should work similarly as they do on scRNA-seq data.
If you are running in to any issues or would like a discussion about use of doublet detecting softwares on snRNA-seq data, please feel free to reach out.



.. _snatac:

snATAC-seq
----------
Demultiplexing Softwares
^^^^^^^^^^^^^^^^^^^^^^^^
Demultiplexing snATAC-seq data can be done with the current demultiplexing softwares. 
However, we note that it is much more memory and time consumptive than scRNA-seq.
Additionally, the SNPs should be filtered by SNPs overlapping peak locations instead of exon or gene locations.
You may even want to filter the SNPs further if you still have many after filtering on minor allele frequency and peak location.
We typically aim for ~250,000 SNP.
Regardless, since UMI tags aren't used for snATAC-seq data, demultiplexing can take a lot of memory and time.

In addtion, the following flags are required for each of the following softwares to effectively process snATAC-seq data.

**Souporcell**

- ``--no_umi True``

Doublet Detecting Softwares
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Technically, doublet detecting softwares cannot be applied to snATAC-seq data as they rely on the unique transcriptomes of each cell type to identify heterotypic doublets.
However, if snATAC-seq peaks can be made into a scRNA-seq-like matrices (`i.e.` by linking peaks to genes or some other method), the doublet detecting softwares outlined in this workflow could be applied to snATAC-seq data.
This has been shown previously by `SnapATAC <https://www.nature.com/articles/s41467-021-21583-9>`__ and `ArchR <https://www.nature.com/articles/s41467-021-21583-9>`__ has a method built in that uses a very similar method to :ref:`Scrublet <Scrublet-docs>`.
`AMULET <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02469-x>`__ is another doublet detecting method that has been developed specifically for snATAC-seq data.



Combined snRNA-seq + snATAC-seq
-------------------------------
.. _combo-demultiplex:

Demultiplexing Softwares
^^^^^^^^^^^^^^^^^^^^^^^^
We have noticed a higher percentage of ambient RNA from our combined snRNA-seq + scATAC-seq experiments as compared to our scRNA-seq (we haven't tested multiplexed snRNA-seq in our hands) but similar snATAC-seq ambient DNA estimations as detected with :ref:`Souporcell<Souporcell-docs>`
Therefore, we recommend running :ref:`Souporcell<Souporcell-docs>`, if only to estimate the ambient RNA in your multiplexed pool.
If you are doing the demultiplexing with the snRNA-seq results, please see the :ref:`snRNA-seq Section <snrna>`.
If you are using the snATAC-seq data for demultiplexing, please see the :ref:`snATAC-seq Section <snatac>`.


Doublet Detecting Softwares
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Doublet detecting softwares for the combined snRNA-seq + snATAC-seq should work similarly to the doublet detecting softwares for each assay separately (snRNA-seq and snATAC-seq).
However, as noted int the :ref:`Demultiplexing Softwares Section <combo-demultiplex>` above, we have observed much higher ambient RNA percentages than for either assay run separately.
..Our results (**CITATION**) indicate that increased ambient RNA showed a slight decrease in the MCC and balanced accuracy. However, we did not simulate up to the level of ambient RNA percent that we have observed using this assay.