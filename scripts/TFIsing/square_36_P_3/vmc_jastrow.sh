#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham ising --boundary peri --ham_dim 2 --L 6 --h 3 --net jas --seed 123 --show_progress
