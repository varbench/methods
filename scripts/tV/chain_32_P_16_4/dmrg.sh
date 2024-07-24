#!/bin/sh

cd programs/dmrg_itensors || exit
julia --project tv_1d.jl --peri --L 32 --V 4 --Nf 16 --max_B 200 --seed 123
