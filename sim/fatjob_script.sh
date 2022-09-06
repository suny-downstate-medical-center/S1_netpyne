#!/bin/bash
#SBATCH --nodes=8            # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=12:00:00               # time limits: 1 hour
#SBATCH --partition=g100_usr_bmem
#SBATCH --account=icei_H_King

srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch1/v12_batch1_0_cfg.json netParams=../data/v12_batch1/v12_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch2/v12_batch2_0_cfg.json netParams=../data/v12_batch2/v12_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch3/v12_batch3_0_cfg.json netParams=../data/v12_batch3/v12_batch3_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch64/v12_batch64_0_cfg.json netParams=../data/v12_batch64/v12_batch64_netParams.py



