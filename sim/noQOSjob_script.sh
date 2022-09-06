#!/bin/bash
#SBATCH --nodes=8            # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=24:00:00               # time limits: 1 hour
#SBATCH --partition=g100_usr_prod
#SBATCH --account=icei_H_King

srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch1/v12_batch1_0_cfg.json netParams=../data/v12_batch1/v12_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch2/v12_batch2_0_cfg.json netParams=../data/v12_batch2/v12_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch3/v12_batch3_0_cfg.json netParams=../data/v12_batch3/v12_batch3_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch4/v12_batch4_0_cfg.json netParams=../data/v12_batch4/v12_batch4_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch5/v12_batch5_0_cfg.json netParams=../data/v12_batch5/v12_batch5_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch6/v12_batch6_0_cfg.json netParams=../data/v12_batch6/v12_batch6_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch7/v12_batch7_0_cfg.json netParams=../data/v12_batch7/v12_batch7_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch8/v12_batch8_0_cfg.json netParams=../data/v12_batch8/v12_batch8_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch9/v12_batch9_0_cfg.json netParams=../data/v12_batch9/v12_batch9_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch10/v12_batch10_0_cfg.json netParams=../data/v12_batch10/v12_batch10_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch11/v12_batch11_0_cfg.json netParams=../data/v12_batch11/v12_batch11_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch12/v12_batch12_0_cfg.json netParams=../data/v12_batch12/v12_batch12_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch13/v12_batch13_0_cfg.json netParams=../data/v12_batch13/v12_batch13_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch14/v12_batch14_0_cfg.json netParams=../data/v12_batch14/v12_batch14_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch15/v12_batch15_0_cfg.json netParams=../data/v12_batch15/v12_batch15_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v12_batch16/v12_batch16_0_cfg.json netParams=../data/v12_batch16/v12_batch16_netParams.py

