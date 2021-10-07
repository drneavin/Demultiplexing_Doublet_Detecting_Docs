.. Demultiplexing & Doublet Detection Workflow documentation master file, created by
   sphinx-quickstart on Sun May 30 20:33:27 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.




Demultiplexing & Doublet Detection Workflow
=======================================================================
.. sidebar:: Number and Percent Doublets per Number Droplets Captured

    .. image:: https://user-images.githubusercontent.com/44268007/120289463-78db8200-c300-11eb-8317-9614efff9838.png


As droplet-based single cell technologies have advanced, increasingly larger sample numbers have been used to answer research questions at single cell resolution.
With a larger number of droplets captured, there is an increase in the proportion of the droplets that are doublets (*Figure 1*).

This has been made possible because, as the droplet-based capture technologies have been optimized, methods to pool and then demultiplex samples - assign droplets to each individual in the pool - have been developed.
These multiplexing methods clearly decrease cost and time of scRNA-seq experiments. 
If left in the dataset, doublets can significantly impact scientific conclusions such identifying spurious cell trajectories or false novel cell types. 
Therefore, it's crucial to effectively clean datasets prior to downstream analyses.

In addition to demultiplexing softwares, there are also doublet detecting softwares that use the transcriptional profiles of droplets to identify doublets by simulating doublets.



.. toctree::
   :maxdepth: 2
   :caption: General
   :hidden:
   
   Background
   Installation
   DataPrep
   SoftwareSelections
   Doublet Estimation Calculator <test.html#http://>
   OtherSingleCellData


.. toctree::
   :maxdepth: 2
   :caption: Demultiplexing Tutorials
   :hidden:

   DemultiplexingSoftwares
   Demuxlet
   Freemuxlet
   scSplit
   Souporcell
   Vireo


.. toctree::
   :maxdepth: 3
   :caption: Doublet Detecting Tutorials
   :hidden:

   DoubletDetectingSoftwares
   DoubletDecon
   DoubletDetection
   DoubletFinder
   Scds
   scDblFinder
   Scrublet
   Solo


- `Doublet Estimation Calculator`__ for testing

__ test.html

.. Reference
.. ==================



Support
==================
If you're having trouble with any part of the Demultiplexing and Doublet Detecting Pipeline, feel free to submit an `issue <https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/issues>`_.
