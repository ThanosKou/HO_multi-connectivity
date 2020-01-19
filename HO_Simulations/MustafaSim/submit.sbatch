#!/bin/bash

#SBATCH --mem=2GB
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%a_%A.out
#SBATCH --error=logs/%a_%A.err
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=apk337@nyu.edu


module purge
module load matlab/2019a

cd /scratch/$USER/HO_multi-connectivity/HO_Simulations/MustafaSim


echo "$SLURM_ARRAY_TASK_ID"


matlab -nodisplay -r "SimulationLOS, exit"


