VMC and ED implemented using [NetKet](https://github.com/netket/netket).

The code requires Python 3.10.4 . `requirements.txt` contains a simple set of dependencies to reproduce the results. It does not include CUDA, and it is enough for quick experiments on small lattices. Create a virtual environment and use `pip install -r requirements.txt` to install them.

`vmc.py` runs the VMC with simple ansatzes such as Jastrow, RBM, and RNN, which are used as baselines in the VarBench dataset.

`ed.py` runs the naive ED using NetKet. It does not utilize the symmetries, so we only use it on small lattices with <= 20 sites.

`args_parser.py` contains all the configurations.

Alternatively, `conda_env.yaml` contains the dependencies including CUDA, which is usually needed for VMC on larger lattices. It also includes [lattice-symmetries](https://github.com/twesterhout/lattice-symmetries) to implement ED on larger lattices. Use `conda env create -f conda_env.yaml` to install them.

`ed_ls.py` runs the ED using lattice-symmetries.
