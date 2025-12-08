#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --job-name=movies
#SBATCH --partition=short
#SBATCH -M carrington
##SBATCH -A spacephysics
#SBATCH --nodes=1
#SBATCH -c 16                  # CPU cores per task
#SBATCH -n 1                  # number of tasks (4xnodes)
#SBATCH --mem=32G # memory per node


#for i in 1 2 3 4 5 6 7 8 9 10
for i in 7
do

/proj/group/spacephysics/ffmpeg/ffmpeg -threads 16 -y -f image2 -r 10 -start_number 1073 -i fluxrope_${i}.0/FHA_multimap_y_v_${i}.0_%07d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -r 25 -f mp4 -vcodec libx264 -pix_fmt yuv420p FHA_FXO_cusp_${i}.mp4
#/proj/group/spacephysics/ffmpeg/ffmpeg -threads 16 -y -f image2 -r 10 -start_number 1056 -i fluxrope_${i}.0/FHA_multimap_y_v_${i}.0_%07d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,transpose=2" -r 25 -f mp4 -vcodec libx264 -pix_fmt yuv420p FHA_FXO_${i}_cusp_vertical.mp4

done
