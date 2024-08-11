import numpy as np
import utils
import sys
import config
import lattice_symmetries as ls
import scipy.linalg
import observables
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--j2", type=float, default=0.0)
parser.add_argument("--lattice", type=str, choices=["square4x4", "triangle4x4", "kagome2x3", "square6x4"])
parser.add_argument("--log2samples", type=int, default=36)
parser.add_argument('--symmetry', type=int, default=0)
args = parser.parse_args()


opt_config = config.opt_parameters(args)
H = opt_config.hamiltonian(**opt_config.ham_params_dict)

circuit = opt_config.circuit(**opt_config.circuit_params_dict)

projector = opt_config.projector(H.bonds, **opt_config.proj_params_dict)


obs = None
if not opt_config.with_mpi or rank == 0:
    obs = observables.Observables(opt_config, H, circuit, projector)

opt = opt_config.optimizer(H, circuit, projector, obs, opt_config.algorithm, opt_config, opt_config.opt_params_dict)


print(opt.optimize())
