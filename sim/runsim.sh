#!/bin/bash
#SBATCH --account=def-bojana
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16    # There are 24 CPU cores on Cedar GPU nodes
#SBATCH --mem=500GB  #500GB  
#SBATCH --mail-user=monica.bellvila@mail.utoronto.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE

module load neuron/7.8.2 mpi4py

srun ./x86_64/special -mpi -python init.py
