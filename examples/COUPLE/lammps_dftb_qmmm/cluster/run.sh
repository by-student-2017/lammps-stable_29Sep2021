#!/bin/bash

# Usage
# bash ./run.sh

NCPUs=1
NGPUs=1 # for DFTB+

lammps_adress=$HONE/lammps-stable_29Sep2021 # Usually, Linux setting (e.g., ubuntu)
lammps_adress=/mnt/d/lammps-stable_29Sep2021 # D drive on WSL (windows11)

export OMP_NUM_THREADS=${NCPUs}
export MAGMA_NUM_GPUS=${NGPUs}

mpirun -np ${NCPUs} ${lammps_adress}/src/lmp_mpi -v mode file < in.client.lmp & python2 ./../dftb_wrap.py file geometry_tmp.gen &