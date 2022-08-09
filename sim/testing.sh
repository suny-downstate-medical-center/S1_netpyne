#!/bin/bash

mpiexec -n 2 nrniv -python -mpi init.py simConfig=../data/v11_batch1/v11_batch1_0_cfg.json netParams=../data/v11_batch1/v11_batch1_netParams.py
mpiexec -n 2 nrniv -python -mpi init.py simConfig=../data/v11_batch1/v11_batch1_1_cfg.json netParams=../data/v11_batch1/v11_batch1_netParams.py
