Installation & Preparation
==========================
Installation should be pretty painless (we hope).
Just download the singluarity image with::

	wget ...


Test Dataset
------------
We have shipped the image with a **small** dataset that can be used to test the softwares of interest.

.. WARNING:: 
	The test dataset may not produce real-world results due to small size - especially for doublet detecting softwares

You can get the dataset from singularity image by running::

	singularity exec --bind <bind_path> ...

where ``<bind_path>`` is a directory somewhere above your current directory
