.. |br| raw:: html

   <br />

Overview of Doublet Detecting Softwares
===========================================

Transcrition-based doublet detection softwares use the transcriptomic profiles in each cell to predict whether that cell is a singlet or doublet.
Most methods simulate doublets by adding the transcritional profiles of two droplets in your pool together.
Therefore, these approaches assume that only a small percentage of the droplets in your dataset are doublets.
The table bellow provides a comparison of the different methods.

+--------------------------------------------------+------------------------------------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| Doublet Detecting Software                       | .. centered:: QC Filtering Required      | .. centered:: Requires Pre-clustering    | Doublet Detecting Method                                                                                                                                  |
+==================================================+==========================================+==========================================+===========================================================================================================================================================+
| :ref:`DoubletDecon <DoubletDecon-docs>`          | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_check_mark:|       | Deconvolution based on clusters provided.                                                                                                                 |
+--------------------------------------------------+------------------------------------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`DoubletDetection <DoubletDetection-docs>`  | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_multiplication_x:| | Iterative boost classifier to classify doublets.                                                                                                          |
+--------------------------------------------------+------------------------------------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`DoubletFinder <DoubletFinder-docs>`        | .. centered:: |:heavy_check_mark:|       | .. centered:: |:heavy_multiplication_x:| | Identify ideal cluster size and call expected number of droplets with highest number of simulated doublet neighbors as doublets.                          |
+--------------------------------------------------+------------------------------------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`scDblFinder <scDblFinder-docs>`            | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_multiplication_x:| | Gradient boosted trees trained with number neighboring doublets and QC metrics to classify doublets                                                       |
+--------------------------------------------------+------------------------------------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`Scds <Scds-docs>`                          | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_multiplication_x:| | **cxds**: Uses genes pairs that are typically not expressed in the same droplet to rank droplets based on coexpression of all pairs. |br|                 |
|                                                  |                                          |                                          | **bcds**: Uses highly variable genes and simulated doublets to train a binary classificaiton algorithm and return probability of droplet being a doublet. |
+--------------------------------------------------+------------------------------------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`Scrublet <Scrublet-docs>`                  | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_multiplication_x:| | Identifies the number of neighboring simulated doublets for each droplet and uses bimodal distribution of scores to classify singlets and doublets.       |
+--------------------------------------------------+------------------------------------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`Solo <Solo-docs>`                          | .. centered:: |:heavy_multiplication_x:| | .. centered:: |:heavy_multiplication_x:| | Simulates doublets and fits a two-layer neural network.                                                                                                   |
+--------------------------------------------------+------------------------------------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+

If you don't know which demultiplexing software(s) to run, take a look at our :ref:`Software Selection Recommendations <SoftwareSelection-docs>` based on your dataset or use our **add widget link here**
