#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker9
#SBATCH --output="CAPS_in_w4_BA_COMPUTE_DM_FINAL_TV_SC_angle0.out"

srun matlab -nosplash -nodesktop -nodisplay -r "BA_COMPUTE_DM_FINAL_TV_SC;  exit"
