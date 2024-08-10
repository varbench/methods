import numpy as np
import utils
import sys
import config as cv_module
import lattice_symmetries as ls
import scipy.linalg
import observables


config_file = utils.import_config(sys.argv[1])
config_import = config_file.opt_parameters()

opt_config = cv_module.opt_parameters()
opt_config.__dict__ = config_import.__dict__.copy()

H = opt_config.hamiltonian(**opt_config.ham_params_dict)

circuit = opt_config.circuit(**opt_config.circuit_params_dict)

projector = opt_config.projector(H.bonds, **opt_config.proj_params_dict)


obs = None
if not opt_config.with_mpi or rank == 0:
    obs = observables.Observables(opt_config, H, circuit, projector)

opt = opt_config.optimizer(H, circuit, projector, obs, opt_config.algorithm, opt_config, opt_config.opt_params_dict)


print(opt.optimize())
