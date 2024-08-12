#!/bin/sh

cd programs/dmrg_itensors_c++ || exit
g++ dmrg_triangular_heisenberg_16x16.cc -o dmrg_triangular_heisenberg_16x16
./dmrg_triangular_heisenberg_16x16
