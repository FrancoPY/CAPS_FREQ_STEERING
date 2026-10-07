#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker8
#SBATCH --output="CAPS_in_w5_generateDensityMap_param2ym4.out"

srun matlab -nosplash -nodesktop -nodisplay -r "generateDensityMap;  exit"
