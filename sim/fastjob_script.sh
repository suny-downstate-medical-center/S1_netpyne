#!/bin/bash
#SBATCH --nodes=2            # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=2:00:00               # time limits: 1 hour
#SBATCH --partition=g100_usr_prod
#SBATCH --qos=g100_qos_dbg
#SBATCH --account=icei_H_King
#SBATCH --job-name=v8_batch11_0_0
#SBATCH -o ../data/v8_batch11/v8_batch11_0_0.run
#SBATCH -e ../data/v8_batch11/v8_batch11_0_0.err

srun ./x86_64/special -mpi -python init.py
