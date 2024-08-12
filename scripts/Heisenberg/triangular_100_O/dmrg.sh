#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project heisenberg_2d.jl --L 10 --J2 1 --zero_mag --max_B 4096 --seed 123

cd programs/dmrg_itensors_c++ || exit
g++ dmrg_triangular_heisenberg_10x10.cc -o dmrg_triangular_heisenberg_10x10
./dmrg_triangular_heisenberg_10x10
