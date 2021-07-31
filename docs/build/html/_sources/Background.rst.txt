Background
==========================
While we described some of the best methods in our manuscript for demultiplexing and doublet detecting, we acknowledge that each dataset is different and may have unique characteristics that make other softwares more suited.
Therefore, this workflow is set up to enable the user to choose and run the demultiplexing and doublet detecting analyses of their choice smoothly and efficiently.
We have built all of the software you will need for any of these softwares into a singularity image that can be easily run on most HPCs.
This means that you do not need to install each software separately and provides standardization of softwares across studies and/or collaborations.
We have also built scripts that will easly summarize the results from each software for you - making the assessment of the success of a software faster and easier.
Finally, we provide a simple command that will easily combine the results from each of the individual softwares into a common dataframe and provide summary statistics about that combination.

We try our best to include all the possible methods for demultiplexing and doublet detecting in this image and maintain them up-to-date. 
If you notice a demultiplexing or doublet detecting software for scRNA-seq data that we have not included, please reach out to us.
