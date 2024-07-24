#!/bin/sh

cd programs/vmc_netket || exit
python3 ed.py --ham hubb --boundary peri --ham_dim 1 --L 14 --U 7.74263683 --Nf 4
