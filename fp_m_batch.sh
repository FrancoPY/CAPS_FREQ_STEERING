#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=gamerpcs
#SBATCH --nodelist=worker3
#SBATCH --output="CAPS_in_w1_densitymapInc.out"

srun matlab -nosplash -nodesktop -nodisplay -r "generateDensityMap_Inc;  exit"
