#!/bin/bash

module purge                 ## Purge modules for a clean start
module load Python/3.7.2-GCCcore-8.2.0
module load matplotlib
export PTNOLATEX=1 PTBACKEND=Agg PTNONINTERACTIVE=1 PYTHONPATH=$PYTHONPATH:/proj/ykempf/analysator/:/proj/ykempf/.local PYTHONUSERBASE=/proj/ykempf/.local

python plot_multi_with_xo.py bulk1 1612 1612 7
