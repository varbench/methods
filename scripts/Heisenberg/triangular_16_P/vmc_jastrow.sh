#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham heis_tri --boundary peri --sign mars --ham_dim 2 --L 4 --zero_mag --net jas --seed 123 --show_progress
