#!/bin/bash
#SBATCH --nodes=8            # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=24:00:00               # time limit
#SBATCH --partition=g100_usr_bmem
#SBATCH --account=icei_H_King

srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch1/v13_batch1_0_cfg.json netParams=../data/v13_batch1/v13_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch2/v13_batch2_0_cfg.json netParams=../data/v13_batch2/v13_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch3/v13_batch3_0_cfg.json netParams=../data/v13_batch3/v13_batch3_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch4/v13_batch4_0_cfg.json netParams=../data/v13_batch4/v13_batch4_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch5/v13_batch5_0_cfg.json netParams=../data/v13_batch5/v13_batch5_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch6/v13_batch6_0_cfg.json netParams=../data/v13_batch6/v13_batch6_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch7/v13_batch7_0_cfg.json netParams=../data/v13_batch7/v13_batch7_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch8/v13_batch8_0_cfg.json netParams=../data/v13_batch8/v13_batch8_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch9/v13_batch9_0_cfg.json netParams=../data/v13_batch9/v13_batch9_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch10/v13_batch10_0_cfg.json netParams=../data/v13_batch10/v13_batch10_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch11/v13_batch11_0_cfg.json netParams=../data/v13_batch11/v13_batch11_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch12/v13_batch12_0_cfg.json netParams=../data/v13_batch12/v13_batch12_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch13/v13_batch13_0_cfg.json netParams=../data/v13_batch13/v13_batch13_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch14/v13_batch14_0_cfg.json netParams=../data/v13_batch14/v13_batch14_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch15/v13_batch15_0_cfg.json netParams=../data/v13_batch15/v13_batch15_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v13_batch16/v13_batch16_0_cfg.json netParams=../data/v13_batch16/v13_batch16_netParams.py

