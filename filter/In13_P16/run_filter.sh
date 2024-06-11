#!/bin/bash 

# Define the necessary paths to executables and 
# pseudopotential files
cd $PBS_O_WORKDIR
cwd=$(pwd)
echo $cwd
pseudoDir="$cwd"
executableBaseDir="/home/your_path_to_repo/FilterBSE"
filterXDir="$executableBaseDir/filter"
scratchDir="/scratch/your_scratch/$PBS_JOBNAME"

# OpenMP settings
export OMP_NUM_THREADS=1
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
bash inputFilter.sh $OMP_NUM_THREADS

# Copy all required input and executable files into the scratch directories
cp $cwd/filter_input.par $sFilterDir/input.par
cp $cwd/filter_input.par $hFilterDir/input.par
rm $cwd/*_input.par
# Pseudopotentials
cp $pseudoDir/pot* $sFilterDir/

# Neighbor list and conf
cp $cwd/allNeigh* $sFilterDir/
cp $cwd/conf.par $sFilterDir/

ln -s $filterXDir/Filter.x $sFilterDir/

# Run the filter executable on a single node with parallel openMP threads
cd $sFilterDir
./Filter.x 2> error.dat > run.dat 
#if [$? -ne 0]; then echo "Job failed!"; fi
wait

# Move filter results to home directory and scratch BSE directory 
cp strain* err* *cube run.dat 'eval'.dat conf.dat input.par $hFilterDir/

