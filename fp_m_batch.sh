#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker7
#SBATCH --output="CAPS_in_w1_generatedensitymapBgnd.out"

srun matlab -nosplash -nodesktop -nodisplay -r "generateDensityMap_Bgnd;  exit"
