#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker11
#SBATCH --output="CAPS_in_w7_simulateBgnd6_CAPS_steering_freq_m2.out"

srun matlab -nosplash -nodesktop -nodisplay -r "simulateBgnd6_CAPS_steering_freq;  exit"
