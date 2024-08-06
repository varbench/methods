#!/bin/sh

cd programs/tV_model_bf || exit
python3 run_nqs.py --L 4 --Nf 5 --V 0.01 --j 1 --bf 1 --depth 2 --feat 4
