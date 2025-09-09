#!/bin/bash

export PTNOLATEX=1 PTBACKEND=Agg PTNONINTERACTIVE=1 PYTHONPATH=$PYTHONPATH:/proj/ykempf/analysator_ykempf/:/proj/ykempf/.local PYTHONUSERBASE=/proj/ykempf/.local
module load matplotlib

# Figure at t=1570
python plot_shrimp.py 1560 1613 1 1570

# Stack for animation
python plot_shrimp.py 1560 1613 1

/proj/group/spacephysics/ffmpeg/ffmpeg -threads 16 -y -f image2 -r 5 -start_number 1560 -i FHA_VSC_shrimp_slices_LMN_%07d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -r 25 -f mp4 -vcodec libx264 -pix_fmt yuv420p FHA_VSC_shrimp_slices_LMN.mp4

