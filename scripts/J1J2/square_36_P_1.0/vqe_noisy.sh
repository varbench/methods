#!/bin/sh

cd programs/VQE || exit
python3 code/main.py --log2samples=14 --j2=1.0 --lattice=square4x4 --symmetry=0 && python3 code/main.py --j2=1.0 --lattice=square4x4 --symmetry=1 --log2samples=14
