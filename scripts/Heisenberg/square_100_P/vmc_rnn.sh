#!/bin/sh

cd programs/vmc_netket || exit
python3 vmc.py --ham heis --boundary peri --sign mars --ham_dim 2 --L 10 --net rnn_lstm --layers 2 --features 16 --seed 123 --optimizer adam --max_step 100000 --show_progress
