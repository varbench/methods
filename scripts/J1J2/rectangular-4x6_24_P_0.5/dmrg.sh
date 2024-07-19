#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project heisenberg_2d.jl --peri --L 4 --L2 6 --J2 0.5 --J22 0.5 --zero_mag --max_B 4096 --seed 123
