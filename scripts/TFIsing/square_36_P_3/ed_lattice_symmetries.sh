#!/bin/sh

cd programs/vmc_netket || exit
python3 ed_ls.py --ham ising --boundary peri --ham_dim 2 --L 6 --h 3
