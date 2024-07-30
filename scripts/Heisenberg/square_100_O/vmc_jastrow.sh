#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham heis --boundary open --sign mars --ham_dim 2 --L 10 --zero_mag --net jas --seed 123 --show_progress
