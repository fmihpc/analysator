#!/usr/bin/env bash
if [[ ! $ARRAY_SIZE ]]; then
    ARRAY_SIZE=14 
fi
# cmd argument $1 is used for the file name of the script to run, log file, jobid 

#if sbatch fails as it gets the error code from python it tries to print the log 
#   first with srun (in case the file has not yet updated on the front end)
#   second if srun fails (in case communication failure) it tries to cat it on the frontend (may be empty if file has not been updated pyproject)

sbatch -W --array=1-$ARRAY_SIZE -o "$1" ./testpackage/$1.sh $2 > jobid_$1 || srun cat $1 || cat $1 || echo "cat failed exit code $?"

#in case we do exit 0 successfully
LOG=$(srun cat $1 || cat $1)
if [[ ! $LOG ]]; then
  echo "::warning::Log file could not be read: srun failed and cat returned and empty logfile. Exit code $?" 
  exit 1
fi
echo "$LOG"
#It is possible that the sbatch command above returns exit 0 if only for example 1 of the array jobs failed but not all, in such a case we check sacct
#   It is also possible that the node never ran it and silently failed which is visible on sacct
#Additionally note that JOBID is saved by the run scripts HOWEVER if a nodes silently fails it may not get passed to the github output
JOBID=$(srun grep -Po '\d+' jobid_$1 || grep -Po '\d+' jobid_$1)
if [[ $? -ne 0 ]] || [[ ! $JOBID ]]; then 
  echo "::error::Jobid could not be found, exiting. Exit code $?" 
  exit 1
fi
SACCT_LOG=$(sacct -j $JOBID -o job,state,node,exit | grep FAILED) 
if [[ $SACCT_LOG ]]; then
  echo "::error::Some job failed on a node, try to take a look at the slurm log."
  echo "$SACCT_LOG" 
  exit 1
fi
