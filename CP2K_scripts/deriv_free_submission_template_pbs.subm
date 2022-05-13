#!/bin/bash
#PBS -l walltime=##HOURS##:0:00,select=1:ncpus=##CORESTOT##:mpiprocs=1:ompthreads=1:mem=186gb
#PBS -A st-gpatey-1
#PBS -N ##JOBNAME##
#PBS -e ##JOBNAME##.stde
#PBS -o ##JOBNAME##.stdo

################################################################################
cd ##CDIR##
echo "Running on host: " `hostname`
echo "Current working directory is: " `pwd`

# Set matlab crash directory
module load gcc/5.4.0
module load matlab/R2019b
export MATLAB_PREFDIR=/home/haydensc/scratch/.matlab/R2019b

# Run job
matlab -r '##JOBNAME##' > ##JOBNAME##.optlog

echo "Job completed at `date`"
exit 0