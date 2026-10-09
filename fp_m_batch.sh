#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker8
#SBATCH --output="CAPS_in_w5_simulateInc9_CAPS_steering_freq_param4.out"

srun matlab -nosplash -nodesktop -nodisplay -r "simulateInc9_CAPS_steering_freq;  exit"
