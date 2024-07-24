#!/bin/sh

cd programs/vmc_netket || exit
python3 ed.py --ham heis --boundary peri --ham_dim 1 --L 20 --zero_mag
