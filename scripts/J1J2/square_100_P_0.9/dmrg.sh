#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project heisenberg_2d.jl --peri --L 10 --J2 0.9 --J22 0.9 --zero_mag --max_B 1024 --seed 123
