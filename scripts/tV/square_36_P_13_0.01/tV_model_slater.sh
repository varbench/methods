#!/bin/sh

cd programs/tV_model_bf || exit
python3 run_nqs.py --L 6 --Nf 13 --V 0.01 --j 1 --bf 0
