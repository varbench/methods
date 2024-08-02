#!/bin/sh

cd programs/RNNWavefunctions/2DHeisenberg || exit

python Run_RNNWF_2DHeis.py --J2 1.0 --L 6 --mag_fixed True --Sz 0 --T0 0.25 --lrthreshold 5e-4 --Nannealing 1000 --Nwarmup 1000 --numunits 300 --numsamples 100
python Run_RNNWF_2DHeis.py --J2 1.0 --L 6 --mag_fixed True --RNN_symmetry nosym --numunits 300 --numsamples 100 --T0 0.25 --lrthreshold 5e-4 --Nwarmup 1000 --Nannealing 10000 --lrthreshold_conv 5e-5 --Nconvergence 25000
python Run_RNNWF_2DHeis.py --J2 1.0 --L 6 --mag_fixed True --RNN_symmetry c2dsym --numunits 300 --numsamples 100 --T0 0.25 --lrthreshold 5e-4 --Nannealing 10000 --lrthreshold_conv 5e-5 --Nconvergence 50000 --lrdecaytime_conv 2000

python Run_RNNWF_2DHeis.py --J2 1.0 --L 8 --lr 1e-5 --RNN_symmetry c2dsym --numunits 300 --Nannealing 0 --Nwarmup 0 --Nconvergence 40000 --numsamples 100 --learning_rate_fixed True --lrthreshold_conv 1e-5 

python Run_RNNWF_2DHeis.py --J2 1.0 --L 10 --lr 1e-5 --RNN_symmetry c2dsym --numunits 300 --Nannealing 0 --Nwarmup 0 --Nconvergence 20000 --numsamples 100 --learning_rate_fixed True --lrthreshold_conv 1e-5 

python Run_RNNWF_2DHeis.py --J2 1.0 --L 12 --lr 1e-5 --RNN_symmetry c2dsym --numunits 300 --Nannealing 0 --Nwarmup 0 --Nconvergence 10000 --numsamples 100 --learning_rate_fixed True --lrthreshold_conv 1e-5 
