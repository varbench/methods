#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project ising_2d.jl --peri --L 6 --h 3 --max_B 1024 --seed 123
