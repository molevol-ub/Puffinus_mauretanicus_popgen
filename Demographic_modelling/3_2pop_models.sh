#!/bin/bash

# Bash script to run 2 population models using dadi

# 1. Activate an environment with python3.9+ version

source activate py39

# 2. If you want to test built-in models with dadi-cli (no_mig; split_mig; asym_mig):

for i in {1..10}
do
dadi-cli InferDM --fs ../../Puffinus.SFS.dadi.fs --model split_mig --lbounds 0.0001 0.0001 0 0 --ubounds 5 5 1 5 --output Puffinus.split_nomig.reduced.demog-params --optimizations 20 --cpu 12 --nomisid --p0 0.05 0.05 0.1 0
done

# 3. If you want to test custom models (see Custom_models directory)

for i in {1..10}
do
dadi-cli InferDM --fs ../../Puffinus.SFS.dadi.fs --model vic_delay_size_prebottle --model-file ./Custom_models/vic_delay_size_prebottle.py --lbounds 0.0001 0.0001 0 0 0.00001 0.0001 0.0001 --ubounds 5 5 1 5 5 5 5--output Puffinus.split_nomig.reduced.demog-params --optimizations 20 --cpu 12 --nomisid --p0 0.05 0.05 0.1 0.1 0.1 0.1 0.1
done

# In all cases you can modify the starting parameters (p0) and lower/upper bounds
