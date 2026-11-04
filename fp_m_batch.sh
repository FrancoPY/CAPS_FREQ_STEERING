#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker14
#SBATCH --output="CAPS_in_w7_simulateInc9_CAPS_steering_freq_m2.out"

srun matlab -nosplash -nodesktop -nodisplay -r "simulateInc9_CAPS_steering_freq;  exit"
