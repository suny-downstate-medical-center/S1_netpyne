#!/bin/bash
#SBATCH --nodes=8            # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=12:00:00               # time limits: 1 hour
#SBATCH --partition=g100_usr_bmem
#SBATCH --account=icei_H_King

srun ./x86_64/special -mpi -python init.py simConfig=../data/v11_batch1/v11_batch1_0_cfg.json netParams=../data/v11_batch1/v11_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v11_batch1/v11_batch1_1_cfg.json netParams=../data/v11_batch1/v11_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v11_batch1/v11_batch1_2_cfg.json netParams=../data/v11_batch1/v11_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v11_batch1/v11_batch1_3_cfg.json netParams=../data/v11_batch1/v11_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v11_batch1/v11_batch1_4_cfg.json netParams=../data/v11_batch1/v11_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v11_batch1/v11_batch1_5_cfg.json netParams=../data/v11_batch1/v11_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v11_batch1/v11_batch1_6_cfg.json netParams=../data/v11_batch1/v11_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v11_batch1/v11_batch1_7_cfg.json netParams=../data/v11_batch1/v11_batch1_netParams.py

