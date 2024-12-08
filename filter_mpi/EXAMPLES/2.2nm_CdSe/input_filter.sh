#!/bin/bash

# Defining the necessary path
cwd=$(pwd)

# Set parameters for the grid
# The inputed number of grid points is just an estimate
# the code will determine the minimum number of gridpoints needed
# to cover the nanocrystal. It is okay to make this number small as
# a guess.
nXGrid="40" 
nYGrid="40"
nZGrid="40"
# Convergence of the algorithm depends on the dGrid value. This value will
# always be used by the code and will not be modified. This might be different than
# your filter version.
dGrid="0.8"

# Set options for the pseudopotentials
interpolatePot="0"
useStrain="0"
crystalStructure="wurzite"
outmostMaterial="CdSe"

# Input parameters specific to the filter-diagonalization technique
# nFilterCycles should be 1-128 to take advantage of 128 Perlmutter processors per node
nFilterCycles='16'
mStatesPerFilter='16'
# These parameters control the energy window of the spectrum that should be filtered. 
# It wastes computational effort to place filters inside the gap, so the vbMax should be
# just above the energy of the HOMO and cbMin should be just below the energy of the LUMO
vbMin='-0.31'; vbMax='-0.24'
cbMin='-0.15'; cbMax='-0.06'
setTargets="1"
# These two numbers control how many energy targets are placed in the VB and CB
# They should add up to mStatesPerFilter -> msVB + msCB = mStatesPerFilter
# Because the valence band has a higher density of states than the conduction band, put
# more energy targets in the VB to improve convergence for less computational time
msVB="10"; msCB="6" 
nCheby='2048'
fermiEnergy="-0.18" # CdSe fermi energy differs depending on what potential is used, but somewhere between -0.16 and -0.18
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
getAllStates="1" # If you only want to print out converged states, make this 0
sigmaECut="0.001" # This is the convergence cutoff for what makes a "true eigenstate"
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

