#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker8
#SBATCH --output="CAPS_in_w3_simuateBgnd.out"

srun matlab -nosplash -nodesktop -nodisplay -r "simuateBgnd6_CAPS_steering_freq;  exit"
