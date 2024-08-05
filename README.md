# VarBench Methods

This repository contains all code to help reproduce the results in the [VarBench](https://github.com/varbench/varbench) dataset.

The folder `programs/` contains collections of programs, and `scripts/` contains scripts and input files to run a program for different Hamiltonians. Each script corresponds to a row in the dataset, located at `scripts/<ham>/<lattice>/<method>.sh`, where `<ham>` is the Hamiltonian type, `<lattice>` is the lattice geometry (see `lattice.md` in the dataset repository), and `<method>` is a brief name for the method.

The script is assumed to run on a common Linux distribution, from the root directory of this repository, and outputs energy, sigma, and energy variance of the data row. It should set up a virtual environment with specific versions of all dependencies, and contain the random seed if needed. Alternatively, you may provide instructions on compiling and setting up the virtual environment in the program folder. The script can be omitted only if the program itself is a self-explanatory script for a specific Hamiltonian.

Note that everything in this repository will be published under the Apache-2.0 license unless otherwise stated, and archived on Zenodo.
