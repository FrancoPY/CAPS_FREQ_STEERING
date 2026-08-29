#!/usr/bin/bash
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=gamerpcs
#SBATCH --nodelist=worker3
#SBATCH --output="CAPS_in_w1.out"

srun matlab -nosplash -nodesktop -nodisplay -r "prueba_cluster;  exit"
