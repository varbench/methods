#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project tv_2d.jl --peri --L 6 --V 0.1 --Nf 13 --max_B 4096 --seed 123
