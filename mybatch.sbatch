#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=apk337@nyu.edu
#SBATCH --output=apk_%j.out
  
module purge
module load matlab/2019a

RUNDIR=$SCRATCH/HO_multi-connectivity/HO_Simulations/
cd $RUNDIR

matlab -nodisplay -nodesktop -r "run SimulationLOS.m"
