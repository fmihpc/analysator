#!/bin/bash -l
###carrington:

#SBATCH -J BS-FHA
#SBATCH --output=slurm-%x.%j.out
#SBATCH -t 23:00:00
#SBATCH -M carrington
#SBATCH --partition=short
#SBATCH --ntasks=1               # tasks per node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=160G # memory per node
#SBATCH --array=0-18%7
     # X-Y%Z means "run with $SLURM_ARRAY_TASK_ID from X to Y, on at most Z nodes at a time"
     # run EGL: 0-17
     # run FHA: 0-18, 7-18 for t>=1000
     # run FIA: 0-13
     ### ^ task number indexed by $SLURM_ARRAY_TASK_ID ^

###SBATCH --no-requeue
###SBATCH --exclude=carrington-2



umask 007
ulimit -c unlimited
cd $SLURM_SUBMIT_DIR

t=8

#--------------------------------------------------------------------
#---------------------DO NOT TOUCH-----------------------------------
nodes=$SLURM_NNODES

## Carrington: has 2x16 cores (i think)
#cores_per_node=32

## Vorna: has 2 x 8 cores
cores_per_node=16

# Hyperthreading
ht=2
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t

#--------------------------------------------------------------------



#export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
#module load Python/3.7.2-GCCcore-8.2.0

NPROC=64
RUN=FHA

module list
time srun python utils/biot_savart.py -nproc $NPROC -task $SLURM_ARRAY_TASK_ID -run $RUN
#time srun sg grp-datacloud-spacephysics -c "python /wrk-vakka/users/horakons/carrington/utils/biot_savart.py" -nproc $NPROC -task $SLURM_ARRAY_TASK_ID -run $RUN             # Jonas suggested using sg to handle weird group permission issue. Instead we changed the permissions on the .vlsv files

echo Job complete.
