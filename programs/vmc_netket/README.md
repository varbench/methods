VMC and ED implemented using [NetKet](https://github.com/netket/netket).

`requirements.txt` contains a simple set of dependencies to reproduce the results. It does not include CUDA, but it is enough for ED using NetKet. Create a virtual environment and use `pip install -r requirements.txt` to install them.

`vmc.py` runs the baseline VMC with simple ansatzes such as Jastrow, RBM, and RNN.

`ed.py` runs the naive ED using NetKet. It does not utilize the symmetries, so we only use it for small lattices with <= 20 sites.

`args_parser.py` contains all the configurations.

TODO: ED with symmetries implemented using [lattice-symmetries](https://github.com/twesterhout/lattice-symmetries)

`conda env create -n vmc_netket -f conda_env.yaml`

For the following larger lattices
```
tV/chain_32_P_16_1
tV/chain_32_P_16_2
tV/chain_32_P_16_4
tV/square_36_P_13_0.01
tV/square_36_P_13_0.1
tV/square_36_P_13_1
tV/square_36_P_13_10
```
