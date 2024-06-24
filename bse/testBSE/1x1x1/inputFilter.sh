#!/bin/bash

# Defining the necessary path
cwd=$(pwd)

# Grid settings
nXGrid="20"
nYGrid="20"
nZGrid="20"
dGrid="0.8"

setTargets="1"
calcPotOverlap="0"
getAllStates="1"
sigmaECut="0.1"
timeHamiltonian="0"
interpolatePot="0"
setSeed="0"; randSeed="123"
printNorm="0"
printCubes="1"; ncubes="2" 
saveCheckpoints="0"

# Input parameters specific to the filter-diagonalization technique
# Filter settings
mStatesPerFilter='6'
nFilterCycles='8'
newtonLength='512'
vbMin='-0.57'; vbMax='-0.37'
cbMin='-0.35'; cbMax='-0.17'
msVB="3"; msCB="3"
SOFlag="1"
NLFlag="1"

nThreads=$1
fermiEnergy="-0.3"

# Print an input.par file for a filter-diagonalization calculation
cat > $cwd/filter_input.par << EOF
nx = $nXGrid
ny = $nYGrid
nz = $nZGrid
dGrid = $dGrid
mStatesPerFilter = $mStatesPerFilter
nFilterCycles = $nFilterCycles
nCheby = $newtonLength
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
printCubes = $printCubes $ncubes
setSeed = $setSeed 
printNorm = $printNorm
EOF
