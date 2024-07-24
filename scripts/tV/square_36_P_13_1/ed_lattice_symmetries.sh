#!/bin/sh

cd programs/vmc_netket || exit
python3 ed_ls.py --ham tv --boundary peri --ham_dim 2 --L 6 --V 1 --Nf 13
