#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker9
#SBATCH --output="CAPS_in_w5_computeDeMuSC_param2.out"

srun matlab -nosplash -nodesktop -nodisplay -r "computeDeMuSC;  exit"
