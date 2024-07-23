#!/bin/sh

cd programs/vmc_netket || exit
python3 ed.py --ham hubb --boundary peri --ham_dim 1 --L 14 --U 2.7825594 --Nf 4
