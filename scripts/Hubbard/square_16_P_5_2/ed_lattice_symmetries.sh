#!/bin/sh

cd programs/vmc_netket || exit
python3 ed_ls.py --ham hubb --boundary peri --ham_dim 2 --L 4 --U 2 --Nf 5
