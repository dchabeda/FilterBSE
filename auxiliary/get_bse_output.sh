cwd=$(pwd)
echo $cwd
scratchDir="${PSCRATCH}/filter_bse/$1"

## declare an array of directories
declare -a dirs=("bse_2_2" "bse_8_8" "bse_16_16" "bse_24_24" "bse_32_32" "bse_40_60")

## now loop through the above array
for i in "${dirs[@]}"
do
   echo "Grabbing output from $i"
   cd $cwd
   # Def scratch directory  
   sBSEDir="$scratchDir/$i"
   # Make necessary home directory    
   hBSEDir="$cwd/$i"
   mkdir -p $hBSEDir

   cd $sBSEDir
   # Move filter results to home directory and scratch BSE directory 
   cp err* run*.dat conf* input.par *eval.par *OS* bs-c* exciton* $hBSEDir/
done

