#!/bin/bash

# Usage
# bash ./run.sh

NCPUs=1

lammps_adress=$HONE/lammps-stable_29Sep2021 # Usually, Linux setting (e.g., ubuntu)
lammps_adress=/mnt/d/lammps-stable_29Sep2021 # D drive on WSL (windows11)

export OMP_NUM_THREADS=${NCPUs}

mpirun -np ${NCPUs} ${lammps_adress}/src/lmp_mpi -in set_in_water.lmp