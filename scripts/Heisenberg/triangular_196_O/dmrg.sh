#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project heisenberg_2d.jl --L 14 --J2 1 --zero_mag --max_B 1024 --seed 123
