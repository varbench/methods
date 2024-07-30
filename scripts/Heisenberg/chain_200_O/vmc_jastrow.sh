#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham heis --boundary open --sign mars --ham_dim 1 --L 200 --zero_mag --net jas --seed 123 --show_progress
