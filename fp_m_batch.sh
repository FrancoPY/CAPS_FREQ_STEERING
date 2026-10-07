#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker9
#SBATCH --output="CAPS_in_w4_beamform_to_bf_COMPROBAR.out"

srun matlab -nosplash -nodesktop -nodisplay -r "beamform_to_bf;  exit"
