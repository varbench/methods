#!/bin/sh

cd programs/tV_model_bf || exit
python3 run_nqs.py --L 6 --Nf 13 --V 10.0 --j 1 --bf 1 --depth 2 --feat 6
