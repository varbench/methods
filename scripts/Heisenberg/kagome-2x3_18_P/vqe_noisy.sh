#!/bin/sh

cd programs/VQE || exit
python3 code/main.py --log2samples=14 --j2=0.0 --lattice=kagome3x2 --symmetry=0 && python3 code/main.py --j2=0.0 --lattice=kagome3x2 --symmetry=1 --log2samples=14
