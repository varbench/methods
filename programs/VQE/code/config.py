import numpy as np
import circuits
import hamiltonians
import optimizers
import observables
from scipy.optimize import differential_evolution, minimize
import projector
import utils
import lattice_symmetries as ls
import sys
import os

class opt_parameters:
    def __init__(self):
        self.with_mpi = False
        j2 = float(sys.argv[2])
        n_trial = int(sys.argv[3])

        self.path_to_logs = './j2_{:.3f}_trial_{:d}_samples_{:d}/'.format(j2, n_trial, int(sys.argv[4]))
        os.makedirs(self.path_to_logs, exist_ok=True)
        self.mode = 'continue'
        
        self.target_norm = 0.98 if j2 > 1 else 0.49
        self.lagrange = False
        self.Z = 300.

        self.test = False
        self.reg = 'diag'
        self.state_target = 0

        ### setting up geometry and parameters ###
        self.Lx, self.Ly, self.subl = 4, 4, 1
        self.su2 = True
        self.BC = 'PBC'
        self.spin = 0
        self.noise = False;
        self.noise_p = 0.0;


        self.basis = ls.SpinBasis(ls.Group([]), number_spins=self.Lx * self.Ly, \
            hamming_weight=(self.Lx * self.Ly) // 2 + self.spin if self.su2 else None)
        self.basis.build()
        
        ### setting up symmetries ###
        if int(sys.argv[5]) == 1:
            self.symmetries = [
               utils.get_Cx_symmetry_map(self.Lx, self.Ly, basis=self.basis, su2=self.su2), \
                utils.get_y_symmetry_map(self.Lx, self.Ly, basis=self.basis, su2=self.su2), \
                utils.get_x_symmetry_map(self.Lx, self.Ly, basis=self.basis, su2=self.su2), \
                utils.get_rot_symmetry_map(self.Lx, self.Ly, basis=self.basis, su2=self.su2)
            ]
            self.eigenvalues = [1, 1, 1, 1]
            self.sectors = [0, 0, 0, 0]
            self.degrees = [2, self.Ly, self.Lx, 4]
        else:
            self.symmetries = []
            self.eigenvalues = []
            self.sectors = []
            self.degrees = []

        self.unitary_no = np.ones((self.Lx * self.Ly, self.Lx * self.Ly))
        self.unitary = self.unitary_no

        self.hamiltonian = hamiltonians.HeisenbergSquare
        self.ham_params_dict = {'n_qubits' : self.Lx * self.Ly, \
                                'su2' : self.su2, \
                                'basis' : self.basis, \
                                'Lx' : self.Lx, 'Ly': self.Ly, \
                                'j_pm' : 1.0, \
                                'j2': j2, \
                                'xBC' : 'PBC', \
                                'yBC' : self.BC, \
                                'symmetries' : [s[0] for s in self.symmetries], \
                                'permutations' : [s[1] for s in self.symmetries], \
                                'sectors' : self.sectors, \
                                'spin' : self.spin, \
                                'unitary' : self.unitary, \
                                'workdir' : self.path_to_logs, \
                                'state_target' : self.state_target
                                }


        self.dimerization = [(2 * i, 2 * i + 1) for i in range(self.Lx * self.Ly // 2)]
        self.circuit = circuits.SU2_symmetrized
        self.circuit_params_dict = {'Lx' : self.Lx, \
                                    'Ly' : self.Ly, \
                                    'subl' : self.subl, \
                                    'spin' : self.spin, \
                                    'basis' : self.basis, \
                                    'config' : self, \
                                    'unitary' : self.unitary, \
                                    'BC' : self.BC}
        

        self.projector = projector.ProjectorFull

        self.proj_params_dict = {'n_qubits' : self.Lx * self.Ly, \
                                 'su2' : self.su2, \
                                 'basis' : self.basis, \
                                 'generators' : self.symmetries, \
                                 'eigenvalues' : self.eigenvalues, \
                                 'degrees' : self.degrees}

        self.observables = [observables.neel_order(self.Lx, self.Ly, self.basis, self.su2), \
                            observables.stripe_order(self.Lx, self.Ly, self.basis, self.su2), \
                           ]


        self.optimizer = optimizers.Optimizer
        self.algorithm = optimizers.natural_gradiend_descend
        self.write_logs = True

        self.opt_params_dict = {'lr' : 1e-3}

        #### stochastic parameters ####
        self.N_samples = 2 ** int(sys.argv[4])
        self.SR_eig_cut = 3e-3
        self.SR_diag_reg = 0.
        self.SR_scheduler = False
        return
