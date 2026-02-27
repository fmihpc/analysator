#!/bin/bash
#Make container.def
if [ -f ./container.def ]; then 
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

        zypper -n --no-gpg-checks install python312 python312-pip python312-setuptools python312-devel git

        update-alternatives --install /usr/bin/python python /usr/bin/python3.12 3
        update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.12 3
        update-alternatives --install /usr/bin/pip pip /usr/bin/pip3.12 3
        update-alternatives --install /usr/bin/pip3 pip3 /usr/bin/pip3.12 3

        export PYTHONNOUSERSITE=1
        pip install git+https://github.com/alhom/analysator-backends/releases/tag/v0.0.2
        pip install --editable git+https://github.com/fmihpc/analysator#egg=analysator"
  EOF

  #Build the container
  module load LUMI systools
  singularity build container.sif container.def
else
  export PYTHONNOUSERSITE=1
  singulairty exec container.def bash
fi



