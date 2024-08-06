# tV_model

Copied from https://github.com/imi-hub/tV_model

This repository contains the source code for the symmetrized backflow neural network ([arXiv:2406.09077](https://arxiv.org/abs/2406.09077)).

The method has been implemented in Python using NetKet package and for ED results QuSpin package. Use the following to install the requirements:
```sh
conda env create -f conda_env.yaml
conda activate tV_model_bf
conda install -c weinbe58 --no-deps quspin=0.3.7=py310h36496fc_0
```
