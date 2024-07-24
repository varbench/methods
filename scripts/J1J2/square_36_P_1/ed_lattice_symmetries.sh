#!/bin/sh

cd programs/vmc_netket || exit
python3 ed_ls.py --ham j1j2 --boundary peri --ham_dim 2 --L 6 --J2 1 --zero_mag
