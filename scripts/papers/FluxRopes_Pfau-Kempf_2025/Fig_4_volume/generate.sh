#!/bin/bash

export PTNOLATEX=1 PTBACKEND=Agg PTNONINTERACTIVE=1 PYTHONPATH=$PYTHONPATH:/proj/ykempf/analysator_ykempf/:/proj/ykempf/.local PYTHONUSERBASE=/proj/ykempf/.local
module load matplotlib

python volume.py 1100 1613 25 2 1.0 10.1 0.5 2
