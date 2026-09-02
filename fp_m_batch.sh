#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=thinkstation
#SBATCH --nodelist=worker7
#SBATCH --output="CAPS_in_w1_simulateBgnd.out"

srun matlab -nosplash -nodesktop -nodisplay -r "simulateBgnd6_CAPS_steering_freq;  exit"
