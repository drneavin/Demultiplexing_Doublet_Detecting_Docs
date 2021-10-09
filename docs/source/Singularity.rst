Notes About Singularity Images
===============================

Singularity images effectively store an operating system with files, softwares etc. that can be easily transported across different operating systems - ensuring reproducibility.
Most HPCs have singularity installed making it easy to implement.
There are some tips and tricks we have identified through using singularity images that we thought might help new users.

Tips and Tricks
---------------
1. Error: File Not Found
^^^^^^^^^^^^^^^^^^^^^^^^
  **Reason**

  Singularity only loads the directories directly downstream from where you execute the singularity command.
  If any of the files that need to be accessed by the command are not downstream of the that location, you will receive an error similar to this one:

  .. code-block:: bash

    Failed to open file "/path/to/readfile.tsv" : No such file or directory

  If you then check for that file:

  .. code-block:: bash

    ll /path/to/readfile.tsv

  We can see that the  file does truly exist:

  .. code-block:: bash

    -rw-rw-r-- 1 user group 70636291 Dec 21  2020 /path/to/readfile.tsv

  **Solution**

  The easiest solution to this problem is to "bind" a path upstream of all the files that will need to be accessed by your command:

  .. code-block:: bash

    singularity exec --bind /path Demuxafy.sif ...


If you don't have access to Singularity on your HPC, you can ask your HPC administrators to install it (see the `Singularity page <https://sylabs.io/guides/3.0/user-guide/quick_start.html>`__)