#!/bin/bash

export PTNOLATEX=1 PTBACKEND=Agg PTNONINTERACTIVE=1 PYTHONPATH=$PYTHONPATH:/proj/ykempf/analysator_ykempf/:/proj/ykempf/.local PYTHONUSERBASE=/proj/ykempf/.local
module load matplotlib

python latitude_mapping.py 1073 1613 1

