#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham ising --boundary peri --ham_dim 1 --L 32 --h 1 --net gcnn --layers 2 --features 32 --seed 123 --show_progress
