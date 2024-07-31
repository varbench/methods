#!/bin/sh

cd programs/SpinlesstV-LCT-INT || exit
make
cd -

cd "$(dirname "$0")"/lct_int_inputs || exit
singularity exec ../../../../programs/SpinlesstV-LCT-INT/main.sif mpirun -np 16 /project/SpinlesstV-LCT-INT-PBC/bin/main -a 10 params.in
