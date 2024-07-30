#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham heis_tri --boundary open --sign mars --ham_dim 2 --L 14 --zero_mag --net jas --seed 123 --show_progress
