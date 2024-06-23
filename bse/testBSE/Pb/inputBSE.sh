#!/bin/bash

# Defining the necessary path
cwd=$(pwd)

timingSpecs="0"
printFPDensity="0"; 
saveCheckpoints="0"
n_FP_density="2"
calcDarkStates="0"


# Input parameters specific to the filter-diagonalization technique
# BSE settings
maxHoleStates='4'
maxElecStates='4'
sigmaECut='0.0001'
epsX="6.1"; epsY="6.1"; epsZ="6.1"
deltaEhole="0.02"
deltaEelec="0.02"


SOFlag="1"
NLFlag="1"

nThreads=$1
fermiEnergy="-0.18"

# Print an input.par file for a filter-diagonalization calculation
cat > $cwd/bse_input.par << EOF
maxHoleStates = $maxHoleStates
maxElecStates =	$maxElecStates
sigmaECut = $sigmaECut
epsX = $epsX
epsY = $epsY
epsZ = $epsZ
nThreads = $nThreads
spinOrbit = $SOFlag
NonLocal = $NLFlag
calcDarkStates = $calcDarkStates
timingSpecs = $timingSpecs
fermiEnergy = $fermiEnergy
deltaEhole = $deltaEhole
deltaEelec = $deltaEelec
EOF
