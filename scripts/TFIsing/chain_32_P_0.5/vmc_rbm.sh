#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham ising --boundary peri --ham_dim 1 --L 32 --h 0.5 --net rbm --seed 123 --show_progress
