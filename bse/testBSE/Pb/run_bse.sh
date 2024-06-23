#!/bin/bash 

# Define the necessary paths to executables and 
# pseudopotential files
cd $PBS_O_WORKDIR
cwd=$(pwd)
echo $cwd
pseudoDir="/home/dchabeda/perovPotsLR/largegap"
executableBaseDir="/home/dchabeda/FilterBSE"
BSEXDir="$executableBaseDir/bse"
scratchDir="/scratch/dchabeda/$PBS_JOBNAME"
filterDir="/scratch/dchabeda/$PBS_JOBNAME/filter"

# OpenMP settings
export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Make necessary home and scratch directories
sBSEDir="$scratchDir/bse"
hBSEDir="$cwd/bse"
if [ -d "$sBSEDir" ]; then rm -Rf $sBSEDir; fi
if [ -d "$hBSEDir" ]; then rm -Rf $hBSEDir; fi
mkdir -p $sBSEDir
mkdir -p $hBSEDir

# Make input.par files
cd $cwd
bash inputBSE.sh $OMP_NUM_THREADS

# Copy all required input and executable files into the scratch directories
cp $cwd/bse_input.par $sBSEDir/input.par
cp $cwd/bse_input.par $hBSEDir/input.par
rm $cwd/*_input.par

# Filter output
cp $filterDir/conf.dat $sBSEDir
ln -s $filterDir/output.dat $sBSEDir

ln -s $BSEXDir/BSE.x $sBSEDir/

# Run the filter executable on a single node with parallel openMP threads
cd $sBSEDir
./BSE.x 2> error.dat > run.dat 
#if [$? -ne 0]; then echo "Job failed!"; fi
wait

# Move filter results to home directory and scratch BSE directory 
cp *cube err* run.dat input.par $hBSEDir/

