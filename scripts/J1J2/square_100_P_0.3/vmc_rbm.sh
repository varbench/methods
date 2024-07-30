#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham j1j2 --boundary peri --sign mars --ham_dim 2 --L 10 --J2 0.3 --zero_mag --net rbm --seed 123 --show_progress
