#!/bin/sh

cd programs/vmc_netket || exit
python3 ed.py --ham tv --boundary peri --ham_dim 2 --L 4 --V 10 --Nf 5
