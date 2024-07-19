#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project heisenberg_2d.jl --peri --L 14 --zero_mag --max_B 512 --seed 123
