#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project tv_2d.jl --peri --L 4 --V 1 --Nf 5 --max_B 4096 --seed 123
