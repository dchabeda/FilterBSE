#!/bin/bash 

# Define the necessary paths to executables and 
# pseudopotential files
cwd=$(pwd)
echo $cwd
executableBaseDir="/global/common/software/m2206/FilterBSE"
pseudoDir="${executableBaseDir}/pots/pots_II_VI/"
filterXDir="$executableBaseDir/filter"
scratchDir="${PSCRATCH}/filter_bse/$1"

# OpenMP settings
export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Make necessary home and scratch directories
sFilterDir="$scratchDir/filter"
hFilterDir="$cwd/filter"
if [ -d "$sFilterDir" ]; then rm -Rf $sFilterDir; fi
if [ -d "$hFilterDir" ]; then rm -Rf $hFilterDir; fi
mkdir -p $sFilterDir
mkdir -p $hFilterDir

# Make input.par files
cd $cwd
bash input_filter.sh $OMP_NUM_THREADS

# Copy all required input and executable files into the scratch directories
cp $cwd/filter_input.par $sFilterDir/input.par
cp $cwd/filter_input.par $hFilterDir/input.par
rm $cwd/*_input.par
# Pseudopotentials
# For single geometry calculation
cp $pseudoDir/pot* $sFilterDir/
#cp $pseudoDir/cubic/SO* $sFilterDir/

# Atomic configuration
cp $cwd/conf.par $sFilterDir/

ln -s $filterXDir/Filter.x $sFilterDir/

# Run the filter executable on a single node with parallel openMP threads
cd $sFilterDir
srun -n 1 -c $OMP_NUM_THREADS --cpu_bind=cores ./Filter.x 2> error.dat > run.dat
#if [$? -ne 0]; then echo "Job failed!"; fi
wait


# Move filter results to home directory and scratch BSE directory 
cp err* run.dat conf.dat input.par eval.dat *.cube $hFilterDir/
