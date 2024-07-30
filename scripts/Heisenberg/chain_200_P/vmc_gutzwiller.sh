#!/bin/sh

../../../programs/vmc_gutzwiller/compile.sh
../../../programs/vmc_gutzwiller/lattice_1d.x <wf.in >wf.out
../../../programs/vmc_gutzwiller/preptj_1d.x <datprep.in >prep.out
../../../programs/vmc_gutzwiller/monte_carlo_tj.x <VMC.in >VMC.out







