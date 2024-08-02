Python codes for the parametrized circuit ansatz.
We employ the hamiltonian variational (HV) and the heuristic RYCNOT ansatz.
The codes do not perform parameter optimization.
Here we simply load the optimal parameters which have been found with the VQE algorithm.
We consider a specific number of blocks (circuit depth) for each case: 24 layers for the HV circuit and 10 layers for the RYCNOT.
These correspond to the best variational energies shown in the main text.

The python code requires qutip package (https://qutip.org), version 4.7.0.
