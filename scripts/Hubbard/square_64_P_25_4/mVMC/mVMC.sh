#!/bin/sh

cd programs/mVMC_RBM || exit
./build.sh
cd -

cd "$(dirname "$0")"
python3 ../../../../programs/mVMC_makedef/make_all.py
cd def_files || exit
mpiexec -np 8 ../../../../../programs/mVMC_RBM/build/src/mVMC/vmc.out -e namelist.def
