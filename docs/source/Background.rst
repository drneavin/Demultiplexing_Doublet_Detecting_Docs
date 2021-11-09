Background
==========================


.. sidebar:: Number and Percent Doublets per Number Droplets Captured

    .. image:: https://user-images.githubusercontent.com/44268007/120289463-78db8200-c300-11eb-8317-9614efff9838.png


As droplet-based single cell technologies have advanced, increasingly larger sample numbers have been used to answer research questions at single cell resolution.
With a larger number of droplets captured, there is an increase in the proportion of the droplets that are doublets (*Figure 1*).

This has been made possible because, as the droplet-based capture technologies have been optimized, methods to pool and then demultiplex samples - assign droplets to each individual in the pool - have been developed.
These multiplexing methods clearly decrease cost and time of scRNA-seq experiments. 
If left in the dataset, doublets can significantly impact scientific conclusions such as identifying spurious cell trajectories or false novel cell types. 
Therefore, it's crucial to effectively clean datasets prior to downstream analyses.

In addition to demultiplexing softwares, there are also doublet detecting softwares that use the transcriptional profiles of droplets to identify doublets by simulating doublets.


.. While we described some of the best methods in our manuscript for demultiplexing and doublet detecting, we acknowledge that each dataset is different and may have unique characteristics that make other softwares more suited.
.. Therefore, this workflow is set up to enable the user to choose and run the demultiplexing and doublet detecting analyses of their choice smoothly and efficiently.
.. We have built all of the software you will need for any of these softwares into a singularity image that can be easily run on most HPCs.
.. This means that you do not need to install each software separately and provides standardization of softwares across studies and/or collaborations.
.. We have also built scripts that will easly summarize the results from each software for you - making the assessment of the success of a software faster and easier.
.. Finally, we provide a simple command that will easily combine the results from each of the individual softwares into a common dataframe and provide summary statistics about that combination.

.. We try our best to include all the possible methods for demultiplexing and doublet detecting in this image and maintain them up-to-date. 
.. If you notice a demultiplexing or doublet detecting software for scRNA-seq data that we have not included, please reach out to us.