Usage in HPC environment
========================


LUMI Container
--------------
Script for building and running a container on LUMI can be found in ``containerization/lumi_container.sh``. 

Simply execute the bash script to build the container, afterwards the same script can be used for launching the container.

Before building the container set the python version for the container export the environment variable ``PY_VERSION``, which defaults to ``PY_VERSION="3.12"``.

For binding folders into the container, set the environment variable ``SINGULARITY_BIND`` before running the container. For example,  ``export SINGULARITY_BIND="/opt,/data:/mnt"``, this will bind /opt to /opt and /data to /mnt in the container.
For more information on singularity usage, see https://docs.sylabs.io/guides/3.0/user-guide/index.html
