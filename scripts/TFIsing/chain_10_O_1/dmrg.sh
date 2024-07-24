#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project ising_1d.jl --L 10 --h 1 --max_B 128 --seed 123
