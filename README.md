# NetPyNE implementation of the somatosensory thalamocortical circuits model

## Description
This code reproduces the simulations for the following paper:

Fernando da Silva Borges,  Joao V.S. Moreira,  Lavinia M. Takarabe,  William W. Lytton,  Salvador Dura-Bernal. **Large-scale biophysically detailed model of somatosensory thalamocortical circuits in NetPyNE**. Frontiers in Neuroinformatics. https://doi.org/10.3389/fninf.2022.884245

We have implemented this highly-detailed and complex model S1 model in NetPyNE, using the data available in the Neocortical Microcircuit Collaboration Portal. We also extended the model by adding thalamic circuits, including 6 distinct thalamic populations with intrathalamic, thalamocortical and corticothalamic connectivity derived from experimental data. Our work provides a widely accessible, data-driven and biophysically-detailed model of the somatosensory thalamocortical circuits that can be utilized as a community tool for researchers to study neural dynamics, function and disease.

### Branches
1. `master`: all data, figures, and codes (6,8 GB)
2. `coreneuron`: only the files needed to run the code (564,7 MB)

## Setup and execution

Requires NEURON with Python and MPI support. 

### NEURON libraries 
1. From /sim run `nrnivmodl mod`. This should create a directory called x86_64. 
2. TIn cfg.py make sure cfg.coreneuron = False
3. To run type: `python batch.py` or `mpiexec -n [num_proc] nrniv -python -mpi init.py`

### CoreNEURON libraries
1. From /sim run `nrnivmodl -coreneuron mod`. This should create a directory called x86_64. 
2. In cfg.py make sure cfg.coreneuron = True
3. To run type: `python batch.py` or `mpirun -n [num_proc] ./x86_64/special -mpi -python init.py`

The resumed code to reproduce the Fig. 7: https://github.com/suny-downstate-medical-center/S1_Thal_NetPyNE_Frontiers_2022.


## Overview of file structure:

* `/sim/init.py`: Main executable; calls functions from other modules. Sets what parameter file to use.

* `/sim/netParams.py`: Network parameters

* `/sim/cfg.py`: Simulation configuration

* `/sim/batch.py`: Run multiple simulations

* `/sim/cells`: source files for the different cell types used in the model; these will be imported into netpyne

* `/sim/mod`: NMODL files containing the ionic channel and synaptic mechanisms used in the model 

* `/data`: where the model and simulation data is stored 

* `/info`: information about the network and all the cells needed to build the microcircuit


For further information please contact: fernandodasilvaborges@gmail.com 
