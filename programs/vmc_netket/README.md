VMC and ED implemented using [NetKet](https://github.com/netket/netket).

`requirements.txt` contains a simple set of dependencies to reproduce the results. It does not include CUDA, but it is enough for ED using NetKet. Create a virtual environment and use `pip install -r requirements.txt` to install them.

`vmc.py` runs the baseline VMC with simple ansatzes such as Jastrow, RBM, and RNN.

`ed.py` runs the naive ED using NetKet. It does not utilize the symmetries, so we only use it for small lattices with <= 20 sites.

`args_parser.py` contains all the configurations.

TODO: ED with symmetries implemented using [lattice-symmetries](https://github.com/twesterhout/lattice-symmetries)

`conda env create -n vmc_netket -f conda_env.yaml`

For the following larger lattices
```
Hubbard/square_16_P_4_10
Hubbard/square_16_P_4_2
Hubbard/square_16_P_4_3.5981
Hubbard/square_16_P_4_4
Hubbard/square_16_P_4_6
Hubbard/square_16_P_4_7.74264
Hubbard/square_16_P_4_8
Hubbard/square_16_P_5_10
Hubbard/square_16_P_5_2.1544
Hubbard/square_16_P_5_2
Hubbard/square_16_P_5_3.5981
Hubbard/square_16_P_5_4
Hubbard/square_16_P_5_6
Hubbard/square_16_P_5_7.74264
Hubbard/square_16_P_5_8
J1J2/rectangular-4x6_24_P_0.5
J1J2/square_36_P_0.3
J1J2/square_36_P_0.4
J1J2/square_36_P_0.5
J1J2/square_36_P_0.6
J1J2/square_36_P_0.7
J1J2/square_36_P_0.8
J1J2/square_36_P_0.9
J1J2/square_36_P_1
TFIsing/square_36_P_3
tV/chain_32_P_16_1
tV/chain_32_P_16_2
tV/chain_32_P_16_4
tV/square_36_P_13_0.01
tV/square_36_P_13_0.1
tV/square_36_P_13_1
tV/square_36_P_13_10
```
