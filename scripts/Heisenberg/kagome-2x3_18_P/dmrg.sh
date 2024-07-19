#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project heisenberg_kagome.jl --L 0 --zero_mag --max_B 512 --seed 123
