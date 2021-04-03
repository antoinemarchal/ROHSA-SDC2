#!/bin/bash
#
#SBATCH --job-name=ROHSA-sky_full_v2
#SBATCH --output=result.txt
#
#SBATCH --ntasks=1 # Run on a single CPU
#SBATCH --time=00:05:00 # Time limit hrs:min:sec 
#SBATCH --mem-per-cpu=400 # Job memory request in Mo
#
#SBATCH --partition=prepost
#SBATCH --output=serial_%j.log # Standard output and error log

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/linkhome/rech/genler01/uzp41ct/ROHSA-SDC2/cfitsio-3.49

date
echo "Running ROHSA-SDC2 on a single CPU core"

./../src/ROHSA-SDC2 ../parameters_sky_full_v2.txt 

date


