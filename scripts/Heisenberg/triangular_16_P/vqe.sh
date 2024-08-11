#!/bin/sh

cd programs/VQE || exit
python3 code/main.py --j2=0.0 --lattice=triangle4x4 --symmetry=0 && python3 code/main.py --j2=0.0 --lattice=triangle4x4 --symmetry=1
