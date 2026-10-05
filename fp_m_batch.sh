#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker8
#SBATCH --output="CAPS_in_w4_simulateBgnd6_CAPS_steering_freq_BA6_BGND11.out"

srun matlab -nosplash -nodesktop -nodisplay -r "simulateBgnd6_CAPS_steering_freq;  exit"
