#!/bin/bash

# Defining the necessary path
cwd=$(pwd)

# Set parameters for the grid
nXGrid="36"
nYGrid="36"
nZGrid="36"
dGrid="0.8"

# Set options for the pseudopotentials
interpolatePot="0"
useStrain="1"
crystalStructure="zincblende"
outmostMaterial="InP"

# Input parameters specific to the filter-diagonalization technique
nFilterCycles='16'
mStatesPerFilter='16'
vbMin='-0.23'; vbMax='-0.19'
cbMin='-0.15'; cbMax='-0.1'
setTargets="1"
msVB="10"; msCB="6"
nCheby='2048'
fermiEnergy="-0.16"
setSeed="0"; randSeed=""

# Options for parallelization
nThreads=$1

# Options for spin-orbit calculation
SOFlag="0"
NLFlag="0"

# Set options for additional output
printNorm="0"
printCubes="1"; 
saveCheckpoints="0"
ncubes="4"
calcPotOverlap="0"
getAllStates="1"
sigmaECut="0.001"
timeHamiltonian="0"

# Print an input.par file for a filter-diagonalization calculation
cat > $cwd/filter_input.par << EOF
nx = $nXGrid
ny = $nYGrid
nz = $nZGrid
dGrid = $dGrid
mStatesPerFilter = $mStatesPerFilter
nFilterCycles = $nFilterCycles
nCheby = $nCheby
VBmin = $vbMin
VBmax = $vbMax
CBmin = $cbMin 
CBmax = $cbMax
nThreads = $nThreads
spinOrbit = $SOFlag
NonLocal = $NLFlag
interpolatePot = $interpolatePot
setTargets = $setTargets $msVB $msCB
calcPotOverlap = $calcPotOverlap
getAllStates = $getAllStates
sigmaECut = $sigmaECut 
timeHamiltonian = $timeHamiltonian
fermiEnergy = $fermiEnergy
setSeed = $setSeed $randSeed
printNorm = $printNorm
printCubes = $printCubes $ncubes
saveCheckpoints = $saveCheckpoints
useStrain = $useStrain
crystalStructure = $crystalStructure
outmostMaterial = $outmostMaterial
EOF
