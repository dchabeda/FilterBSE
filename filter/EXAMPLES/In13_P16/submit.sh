#!/bin/bash

#SBATCH --qos debug
#SBATCH -J InP_test_Filter
#SBATCH -o %j-%N.out
#SBATCH --constraint cpu
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --time 00:30:00
#SBATCH -A m2206
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@berkeley.edu

jobname="InP_test_Filter"

module load PrgEnv-intel/8.5.0
module load cray-fftw/3.3.10.6

bash run_filter.sh $jobname


