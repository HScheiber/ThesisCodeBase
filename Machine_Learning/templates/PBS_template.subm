#!/bin/bash
#PBS -l walltime=##HOURS##:0:00,select=##NODES##:ncpus=##CPUS##:mpiprocs=1:ompthreads=1:mem=##MEMORY##
#PBS -A st-gpatey-1
#PBS -N ##MODEL_NAME##
#PBS -e ##MODEL_NAME##.stde
#PBS -o ##MODEL_NAME##.stdo

# Check on some basics:
echo "Running on host: " `hostname`
echo "Changing to directory from which PBS script was submitted."
cd $PBS_O_WORKDIR
echo "Current working directory is now: " `pwd`


# set EXE environment
module purge
module load Software_Collection/2021
module load gcc
module load openmpi
module load gromacs/2019.6-double
module load matlab/R2021a
mkdir $TMPDIR/.matlab
export MATLAB_PREFDIR=$TMPDIR/.matlab


# Run Job
matlab -nodisplay -r "Bayesian_Optimize_LiX_Parameters('##MODEL_NAME##.inp')" >> ##MODEL_NAME##.log

# Clean up
matlab -r "cleanup_BO_log('##MODEL_NAME##.log')"

exit 0