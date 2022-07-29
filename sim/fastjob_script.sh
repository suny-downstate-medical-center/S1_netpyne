#!/bin/bash
#SBATCH --nodes=2            # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=2:00:00               # time limits: 1 hour
#SBATCH --partition=g100_usr_prod
#SBATCH --qos=g100_qos_dbg
#SBATCH --account=icei_H_King

srun ./x86_64/special -mpi -python init.py simConfig=../data/v11_batch0/v11_batch0_0_cfg.json netParams=../data/v11_batch0/v11_batch0_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v11_batch0/v11_batch0_15_cfg.json netParams=../data/v11_batch0/v11_batch0_netParams.py