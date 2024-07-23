#!/bin/sh

cd programs/vmc_netket || exit
python3 ed.py --ham tv --boundary peri --ham_dim 2 --L 4 --V 1 --Nf 5
