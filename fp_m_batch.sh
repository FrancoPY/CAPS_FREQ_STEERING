#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker9
#SBATCH --output="CAPS_in_w3_correction_beamforming.out"

srun matlab -nosplash -nodesktop -nodisplay -r "beamform_to_bf;  exit"
