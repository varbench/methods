#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham ising --boundary peri --ham_dim 1 --L 10 --h 1 --net jas --seed 123 --show_progress
