#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker7
#SBATCH --output="CAPS_in_w3.1_BA_COMPUTE.out"

srun matlab -nosplash -nodesktop -nodisplay -r "BA_COMPUTE_DM_FINAL;  exit"
