#!/bin/sh

cd programs/vmc_netket || exit
python3 ed.py --ham heis --boundary open --ham_dim 1 --L 20
