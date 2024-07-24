#!/bin/sh

cd programs/vmc_netket || exit
python3 ed_ls.py --ham tv --boundary peri --ham_dim 1 --L 32 --V 2 --Nf 16
