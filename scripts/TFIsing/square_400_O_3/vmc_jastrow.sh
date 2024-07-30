#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham ising --boundary open --ham_dim 2 --L 20 --h 3 --net jas --seed 123 --show_progress
