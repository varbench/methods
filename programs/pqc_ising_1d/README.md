Python codes for the parametrized circuit ansatzes.
We employ the Hamiltonian variational (HV) and the heuristic RYCNOT ansatzes.
The codes do not perform parameter optimization.
Here we simply load the optimal parameters which have been found with the VQE algorithm.
We consider a specific number of blocks (circuit depth) for each ansatz: 24 layers for the HV circuit and 10 layers for the RYCNOT.
These correspond to the best variational energies shown in the dataset.

The code requires Python 3.10.4 and QuTiP 4.7.0 . Use `pip install -r requirements.txt` to install the requirements.
