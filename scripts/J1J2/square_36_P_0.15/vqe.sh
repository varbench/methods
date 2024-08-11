#!/bin/sh

cd programs/VQE || exit
python3 code/main.py --j2=0.15 --lattice=square4x4 --symmetry=0 && python3 code/main.py --j2=0.15 --lattice=square4x4 --symmetry=1
