#!/bin/sh

cd "$(dirname "$0")"
cd ../../../programs/mVMC_RBM/ || exit
./build.sh
cd -
cd mVMC_inputs/ || exit
mpiexec -np 8 ../../../../programs/mVMC_RBM/build/src/mVMC/vmc.out -e namelist.def
