#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project heisenberg_2d.jl --peri --L 6 --J2 0.3 --J22 0.3 --zero_mag --max_B 2048 --seed 123
