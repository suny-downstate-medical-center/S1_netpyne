#!/bin/bash
#SBATCH --job-name=v5_batch12_0_0
#SBATCH -A default
#SBATCH -t 72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH -o ../data/v5_batch12/v5_batch12_0_0.run
#SBATCH -e ../data/v5_batch12/v5_batch12_0_0.err
#SBATCH --mail-user=fernandodasilvaborges@gmail.com
#SBATCH --mail-type=end



source ~/.bashrc
cd /home/ext_fernandodasilvaborges_gmail_/S1_netpyne/sim/
mpirun -n 40 nrniv -python -mpi init.py simConfig=../data/v5_batch12/v5_batch12_0_0_cfg.json netParams=../data/v5_batch12/v5_batch12_netParams.py
wait
                    