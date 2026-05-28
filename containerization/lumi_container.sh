#!/bin/bash
#Make container.def

if [[ "$PY_VERSION" ]]; then
  PY_PACKAGE=$(echo "$PY_VERSION" | sed 's/\.//g' )
else
  PY_VERSION="3.12"
  PY_PACKAGE="312"
fi

if [ ! -f ./container.sif ]; then 
  cat << EOF > container.def
Bootstrap: docker
From: docker.io/opensuse/leap:15.6

%post
   # Continue to install software into the container normally.
   # Building:
       # module load LUMI systools
       # singularity build container.sif container.def
   # Running:
       # For example launching python (one can also launch bash for example)
       # singularity exec container.sif python
       # to bind folders once can use --bind /opt,/data:/mnt
       # this will bind /opt to /opt and /data to /mnt, alternatively one can use environment variables
       # export SINGULARITY_BIND="/opt,/data:/mnt"

        zypper -n --no-gpg-checks install python$PY_PACKAGE python$PY_PACKAGE-pip python$PY_PACKAGE-setuptools python$PY_PACKAGE-devel git

        update-alternatives --install /usr/bin/python python /usr/bin/python$PY_VERSION 3
        update-alternatives --install /usr/bin/python3 python3 /usr/bin/python$PY_VERSION 3
        update-alternatives --install /usr/bin/pip pip /usr/bin/pip$PY_VERSION 3
        update-alternatives --install /usr/bin/pip3 pip3 /usr/bin/pip$PY_VERSION 3

        export PYTHONNOUSERSITE=1
        pip install --extra-index-url https://version.helsinki.fi/api/v4/projects/5244/packages/pypi/simple vlsvrs
        pip install --editable git+"https://github.com/fmihpc/analysator#egg=analysator"
EOF

  #Build the container
  module load LUMI systools
  singularity build container.sif container.def
else
  export PYTHONNOUSERSITE=1
  if [[ ! $SINGULARITY_BIND ]]; then
     echo "Friendly reminder: no folders binded to the container, please use SINGULARITY_BIND environment variable"
  fi 
  singularity exec container.sif bash
fi



