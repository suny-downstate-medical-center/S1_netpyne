#!/bin/bash
#SBATCH --nodes=8            # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=24:00:00               # time limit
#SBATCH --partition=g100_usr_bmem
#SBATCH --account=icei_H_King

srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch1/v14_batch1_0_cfg.json netParams=../data/v14_batch1/v14_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch2/v14_batch2_0_cfg.json netParams=../data/v14_batch2/v14_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch3/v14_batch3_0_cfg.json netParams=../data/v14_batch3/v14_batch3_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch4/v14_batch4_0_cfg.json netParams=../data/v14_batch4/v14_batch4_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch5/v14_batch5_0_cfg.json netParams=../data/v14_batch5/v14_batch5_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch6/v14_batch6_0_cfg.json netParams=../data/v14_batch6/v14_batch6_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch7/v14_batch7_0_cfg.json netParams=../data/v14_batch7/v14_batch7_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch8/v14_batch8_0_cfg.json netParams=../data/v14_batch8/v14_batch8_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch9/v14_batch9_0_cfg.json netParams=../data/v14_batch9/v14_batch9_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch10/v14_batch10_0_cfg.json netParams=../data/v14_batch10/v14_batch10_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch11/v14_batch11_0_cfg.json netParams=../data/v14_batch11/v14_batch11_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch12/v14_batch12_0_cfg.json netParams=../data/v14_batch12/v14_batch12_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch13/v14_batch13_0_cfg.json netParams=../data/v14_batch13/v14_batch13_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch14/v14_batch14_0_cfg.json netParams=../data/v14_batch14/v14_batch14_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch15/v14_batch15_0_cfg.json netParams=../data/v14_batch15/v14_batch15_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v14_batch16/v14_batch16_0_cfg.json netParams=../data/v14_batch16/v14_batch16_netParams.py

