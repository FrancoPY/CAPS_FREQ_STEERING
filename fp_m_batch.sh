#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker9
#SBATCH --output="CAPS_in_w3_correction_computeBA_solo.out"

srun matlab -nosplash -nodesktop -nodisplay -r "BA_COMPUTE_DM_FINAL_SIN_PROMEDIAR;  exit"
