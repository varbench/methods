#!/bin/sh

cd programs/vmc_netket || exit
python3 ed.py --ham ising --boundary open --ham_dim 1 --L 10 --h 1
