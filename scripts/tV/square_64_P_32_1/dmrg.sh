#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project tv_2d.jl --peri --L 8 --V 1 --Nf 32 --max_B 4096 --seed 123
