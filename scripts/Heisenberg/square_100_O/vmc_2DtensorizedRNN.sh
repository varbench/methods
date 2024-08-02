#!/bin/sh

cd programs/RNNWavefunctions/2DHeisenberg || exit

python Run_RNNWF_2DHeis.py --J2 0.0 --L 10 --mag_fixed True --Sz 0 --RNN_symmetry c4vsym --group_character A1 --numunits 200 --Nwarmup 0 --Nannealing 0 --Nconvergence 150000 --lrthreshold_conv 5e-4 --lrdecaytime_conv 5000
