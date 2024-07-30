#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham heis_kag --boundary peri --ham_dim 2 --L 3 --L2 2 --zero_mag --net jas --seed 123 --show_progress
