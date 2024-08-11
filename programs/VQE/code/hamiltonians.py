import os
import os.path
import numpy as np
import scipy as sp
from scipy import sparse
import lattice_symmetries as ls
import utils

sz = np.array([[1, 0], \
               [0, -1]])

sx = np.array(
               [[0, 1], \
               [1, 0]]
              )

sy = np.array([[0, 1.0j], [-1.0j, 0]])

s0 = np.eye(2)
sx_sparse = sp.sparse.csr_matrix(sx)
sy_sparse = sp.sparse.csr_matrix(sy)
sz_sparse = sp.sparse.csr_matrix(sz)
s0_sparse = sp.sparse.csr_matrix(s0)


SS = np.kron(sx, sx) + np.kron(sy, sy) + np.kron(sz, sz)
P_ij = (SS + np.eye(4)) / 2.

SSun = -np.kron(sx, sx) + -np.kron(sy, sy) + np.kron(sz, sz)
P_ijun = (SSun + np.eye(4)) / 2.


class Hamiltonian(object):
    def __init__(self, basis, unitary, n_qubits, su2, symmetries, permutations, sectors, spin, workdir, state_target, **kwargs):
        self.n_qubits = n_qubits
        #self.basis = basis
        self.unitary = unitary
        self.symmetries = symmetries
        self.permutations = permutations
        self.sectors = sectors
        self.spin = spin
   
        ### obtaining ground state in the correct symmetry sector ###
        self.basis = ls.SpinBasis(ls.Group([ls.Symmetry(s, sector=sec) for s, sec in zip(self.symmetries, self.sectors)]), \
                                            number_spins=n_qubits, hamming_weight=n_qubits // 2 + self.spin if su2 else None)#, spin_inversion=-1)
        self.basis.build()
        

        self._matrix, self._terms, self.bonds, self.js = self._get_Hamiltonian_matrix(**kwargs)
        print(self.bonds)


        energy, ground_state = ls.diagonalize(self._matrix, k = 4, dtype=np.complex128)
        self.energies = energy - self.energy_renorm
        state_target = ground_state.T[state_target]
        self.all_states = ground_state
        
        spins = []
        all_bonds = []
        energies_s0 = []
        for i in range(self.n_qubits):
            for j in range(self.n_qubits):
                if i != j:
                    all_bonds.append((i, j))
        total_spin = ls.Operator(self.basis, [ls.Interaction(SS, all_bonds)])        
        for idx, s in enumerate(ground_state.T):
            print(np.dot(s.conj(), total_spin(s)) + 3. * self.n_qubits, n_qubits, 'spin idx')
            spins.append(np.dot(s.conj(), total_spin(s)) + 3. * self.n_qubits)
            if np.abs(spins[-1]) < 1e-2:
                energies_s0.append(self.energies[idx])
        #print(energies_s0[1] - energies_s0[0])
        print(self.energies)
        ### END DEBUG
        for idx, s in enumerate(spins):
            if np.isclose(s / 4, self.spin * (self.spin + 1)):
                break
        print('idx = ', idx)
        
        idx = 0
        ### rewrite ground state in terms of non-symmetrized basis ###
        gs_nonsymm = np.zeros(basis.number_states, dtype=np.complex128)
        state_target_nosymm = np.zeros(basis.number_states, dtype=np.complex128)

        gs_symm = ground_state[:, np.argmin(energy)]  # FIXME
        for i in range(basis.number_states):
            nonsymm_state = basis.states[i]
            rep, character, norm = self.basis.state_info(nonsymm_state)
            if norm != 0.:
                gs_nonsymm[i] = gs_symm[self.basis.index(rep)] * norm * character
                state_target_nosymm[i] = state_target[self.basis.index(rep)] * norm * character


        # gs_nonsymm = gs_nonsymm * np.sqrt(2)  # FIXME FIXME FIXME
        self.ground_state = gs_nonsymm[np.newaxis, :]
        self.state_target = state_target_nosymm

        print(np.vdot(gs_nonsymm, gs_nonsymm))
        #assert np.isclose(np.vdot(gs_nonsymm, gs_nonsymm), 1.0)

        ### finally obtaining the GS in the nonsymmetric basis (provided from config) ###
        self.basis = basis
        self._matrix, self._terms, self.bonds, self.j2s = self._get_Hamiltonian_matrix(**kwargs)   

        #assert np.isclose(np.dot(self._matrix(gs_nonsymm).conj(), gs_nonsymm), np.min(energy))

        #self.ground_state = ground_state.T
        self.nterms = len(self._terms)
        print('ground state energy:', np.min(energy) - self.energy_renorm)
        #print('system gap =', energy[1] - energy[0])

        #print(energy[1] - self.energy_renorm)
        #exit(-1)
        self.gse = np.min(energy)
        return

    def __call__(self, bra, n_term = None, vector=True):
        if n_term is None:
            return self._matrix(bra)
        if vector:
            return self._terms[n_term][0](bra)
        return self._terms[n_term][1]


    def _get_Hamiltonian_matrix(self, **kwargs):
        raise NotImplementedError()


class HeisenbergSquareNNBipartiteOBC(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = -1., j_zz = 1.):
        assert Lx % 2 == 0  # here we only ocnsider bipartite systems 
        assert Ly % 2 == 0

        operator = np.kron(sx, sx) + np.kron(sy, sy) + np.kron(sz, sz)
        n_sites = Lx * Ly

        bonds = []
        for site in range(n_sites):
            x, y = site % Lx, site // Lx

            site_up = ((x + 1) % Lx) + y * Lx
            site_right = x + ((y + 1) % Ly) * Lx

            if x + 1 < Lx:
                bonds.append((site, site_up))
            if y + 1 < Ly:
                bonds.append((site, site_right))
        print('bonds = ', bonds)
        return ls.Operator(self.basis, [ls.Interaction(operator, bonds)]), [ls.Operator(self.basis, [ls.Interaction(operator, [bond])]) for bond in bonds]


class HeisenbergSquareNNBipartitePBC(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = +1., j_zz = 1.):
        assert Lx % 2 == 0  # here we only ocnsider bipartite systems 
        assert Ly % 2 == 0

        operator = j_pm * (np.kron(sx, sx) + np.kron(sy, sy)) + j_zz * np.kron(sz, sz)
        n_sites = Lx * Ly

        bonds = []
        for site in range(n_sites):
            x, y = site % Lx, site // Lx

            site_up = ((x + 1) % Lx) + y * Lx
            site_right = x + ((y + 1) % Ly) * Lx
            bonds.append((site, site_up))
            bonds.append((site, site_right))
        print('bonds = ', bonds)
        return ls.Operator(self.basis, [ls.Interaction(operator, bonds)]), [ls.Operator(self.basis, [ls.Interaction(operator, [bond])]) for bond in bonds]

class HeisenbergSquare(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = +1., j_zz = 1., j2=0., xBC='PBC', yBC = 'PBC'):
        #assert Lx % 2 == 0  # here we only ocnsider bipartite systems
        #assert Ly % 2 == 0

        operator = P_ij
        operator_j2 = P_ij
        operatorun = P_ijun
        operator_j2un = P_ijun

        n_sites = Lx * Ly

        bonds = []
        bonds_j2 = []
        bondsun = []
        bonds_j2un = []

        for site in range(n_sites):
            x, y = site % Lx, site // Lx

            site_up = ((x + 1) % Lx) + y * Lx
            site_right = x + ((y + 1) % Ly) * Lx

            if x + 1 < Lx or xBC == 'PBC':
                if self.unitary[site, site_up] == +1:
                    bonds.append((site, site_up))
                else:
                    bondsun.append((site, site_up))
            if y + 1 < Ly or yBC == 'PBC':
                if self.unitary[site, site_right] == +1:
                    bonds.append((site, site_right))
                else:
                    bondsun.append((site, site_right))

            if not np.isclose(j2, 0.0):
                site_up = ((x + 1) % Lx) + ((y + 1) % Ly) * Lx
                site_right = ((x + 1) % Lx) + ((y - 1) % Ly) * Lx
                if (x + 1 < Lx or xBC == 'PBC') and (y + 1 < Ly or yBC == 'PBC'):
                    if self.unitary[site, site_up] == +1:
                        bonds_j2.append((site, site_up))
                    else:
                        bonds_j2un.append((site, site_up))
                if (x + 1 < Lx or xBC == 'PBC') and (y - 1 >= 0 or yBC == 'PBC'):
                    if self.unitary[site, site_right] == +1:
                        bonds_j2.append((site, site_right))
                    else:
                        bonds_j2un.append((site, site_right))

        print(bonds)
        print(bonds_j2)
        self.energy_renorm = len(bonds) + len(bondsun) + len(bonds_j2) * j2 + len(bonds_j2un) * j2
        return ls.Operator(self.basis, ([ls.Interaction(operator * 2, bonds)] if len(bonds) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2 * 2, bonds_j2)] if len(bonds_j2) > 0 else []) + \
                                       ([ls.Interaction(operatorun * 2, bondsun)] if len(bondsun) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2un * 2, bonds_j2un)] if len(bonds_j2un) > 0 else [])), \
               ([[ls.Operator(self.basis, [ls.Interaction(operator, [bond])]), 2] for bond in bonds] if len(bonds) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2, [bond])]), j2 * 2.] for bond in bonds_j2] if len(bonds_j2) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operatorun, [bond])]), 2] for bond in bondsun] if len(bondsun) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2un, [bond])]), j2 * 2.] for bond in bonds_j2un] if len(bonds_j2un) > 0 else []), \
               bonds + bondsun + bonds_j2 + bonds_j2un, \
               [2] * len(bonds) + [2 * j2] * len(bonds_j2)


class HeisenbergTriangle(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = +1., j_zz = 1., j2=0., xBC='PBC', yBC = 'PBC'):
        #assert Lx % 2 == 0  # here we only ocnsider bipartite systems
        #assert Ly % 2 == 0

        operator = P_ij
        operator_j2 = P_ij
        operatorun = P_ijun
        operator_j2un = P_ijun

        n_sites = Lx * Ly

        bonds = []
        bonds_j2 = []
        bondsun = []
        bonds_j2un = []

        for site in range(n_sites):
            x, y = site % Lx, site // Lx

            site_up = ((x + 1) % Lx) + y * Lx
            site_right = x + ((y + 1) % Ly) * Lx

            if x + 1 < Lx or xBC == 'PBC':
                if self.unitary[site, site_up] == +1:
                    bonds.append((site, site_up))
                else:
                    bondsun.append((site, site_up))
            if y + 1 < Ly or yBC == 'PBC':
                if self.unitary[site, site_right] == +1:
                    bonds.append((site, site_right))
                else:
                    bondsun.append((site, site_right))

            if not np.isclose(j2, 0.0):
                #site_up = ((x + 1) % Lx) + ((y + 1) % Ly) * Lx
                site_right = ((x + 1) % Lx) + ((y - 1) % Ly) * Lx
                #if (x + 1 < Lx or xBC == 'PBC') and (y + 1 < Ly or yBC == 'PBC'):
                #    if self.unitary[site, site_up] == +1:
                #        bonds_j2.append((site, site_up))
                #    else:
                #        bonds_j2un.append((site, site_up))
                if (x + 1 < Lx or xBC == 'PBC') and (y - 1 >= 0 or yBC == 'PBC'):
                    if self.unitary[site, site_right] == +1:
                        bonds_j2.append((site, site_right))
                    else:
                        bonds_j2un.append((site, site_right))

        print(bonds)
        print(bonds_j2)
        self.energy_renorm = len(bonds) + len(bondsun) + len(bonds_j2) * j2 + len(bonds_j2un) * j2
        return ls.Operator(self.basis, ([ls.Interaction(operator * 2, bonds)] if len(bonds) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2 * 2, bonds_j2)] if len(bonds_j2) > 0 else []) + \
                                       ([ls.Interaction(operatorun * 2, bondsun)] if len(bondsun) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2un * 2, bonds_j2un)] if len(bonds_j2un) > 0 else [])), \
               ([[ls.Operator(self.basis, [ls.Interaction(operator, [bond])]), 2] for bond in bonds] if len(bonds) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2, [bond])]), j2 * 2.] for bond in bonds_j2] if len(bonds_j2) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operatorun, [bond])]), 2] for bond in bondsun] if len(bondsun) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2un, [bond])]), j2 * 2.] for bond in bonds_j2un] if len(bonds_j2un) > 0 else []), \
               bonds + bondsun + bonds_j2 + bonds_j2un, \
               [2] * len(bonds) + [2 * j2] * len(bonds_j2)



class HeisenbergKagome12(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = +1., j_zz = 1., j2=0., xBC='PBC', yBC = 'PBC'):
        #assert Lx % 2 == 0  # here we only ocnsider bipartite systems
        #assert Ly % 2 == 0

        operator = P_ij
        operator_j2 = P_ij
        operatorun = P_ijun
        operator_j2un = P_ijun

        n_sites = Lx * Ly

        bonds = [(0, 4), (2, 6), (1, 5), (3, 7), (4, 2), (5, 3), (6, 0), (7, 1), (0, 8), (1, 9), (2, 10), (3, 11), (1, 8), (3, 10), (2, 11), (0, 9), (5, 9), (10, 6), (4, 8), (11, 7), (4, 10), (6, 9), (4, 11), (8, 7)]
        bonds_j2 = [(0, 2), (4, 6), (1, 3), (5, 7), (0, 1), (8, 9), (2, 3), (10, 11), (9, 10), (5, 6), (8, 11), (4, 7)]
        bondsun = []
        bonds_j2un = []

        print(bonds)
        print(bonds_j2)
        self.energy_renorm = len(bonds) + len(bondsun) + len(bonds_j2) * j2 + len(bonds_j2un) * j2
        return ls.Operator(self.basis, ([ls.Interaction(operator * 2, bonds)] if len(bonds) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2 * 2, bonds_j2)] if len(bonds_j2) > 0 else []) + \
                                       ([ls.Interaction(operatorun * 2, bondsun)] if len(bondsun) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2un * 2, bonds_j2un)] if len(bonds_j2un) > 0 else [])), \
               ([[ls.Operator(self.basis, [ls.Interaction(operator, [bond])]), 2] for bond in bonds] if len(bonds) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2, [bond])]), j2 * 2.] for bond in bonds_j2] if len(bonds_j2) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operatorun, [bond])]), 2] for bond in bondsun] if len(bondsun) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2un, [bond])]), j2 * 2.] for bond in bonds_j2un] if len(bonds_j2un) > 0 else []), \
               bonds + bondsun + bonds_j2 + bonds_j2un, \
               [2] * len(bonds) + [2 * j2] * len(bonds_j2)




class HeisenbergHexagon(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = +1., j_zz = 1., j2=0., BC='PBC'):
        operator = P_ij
        operator_j2 = P_ij
        n_sites = 6

        bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5)]
        bonds_j2 = [(0, 2), (2, 4), (4, 0), (1, 3), (3, 5), (1, 5)]
        bondsun = []
        bonds_j2un = []

        self.energy_renorm = len(bonds) + len(bondsun) + len(bonds_j2) * j2 + len(bonds_j2un) * j2
        return ls.Operator(self.basis, ([ls.Interaction(operator * 2, bonds)] if len(bonds) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2 * 2, bonds_j2)] if len(bonds_j2) > 0 else []) + \
                                       ([ls.Interaction(operatorun * 2, bondsun)] if len(bondsun) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2un * 2, bonds_j2un)] if len(bonds_j2un) > 0 else [])), \
               ([[ls.Operator(self.basis, [ls.Interaction(operator, [bond])]), 2] for bond in bonds] if len(bonds) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2, [bond])]), j2 * 2.] for bond in bonds_j2] if len(bonds_j2) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operatorun, [bond])]), 2] for bond in bondsun] if len(bondsun) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2un, [bond])]), j2 * 2.] for bond in bonds_j2un] if len(bonds_j2un) > 0 else []), \
               bonds + bondsun + bonds_j2 + bonds_j2un

class HeisenbergHoneycomb_2x2(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = +1., j_zz = 1., j2=0., BC='PBC'):
        assert ((Lx == 2) and (Ly == 2))
        operator = P_ij
        operator_j2 = P_ij
        n_sites = 8

        bonds = [(0, 1), (0, 5), (0, 7), (1, 4), (1, 6), (2, 3), (2, 5), (2, 7), (3, 4), (3, 6), (4, 5), (6, 7)]
        bonds_j2 = [(0, 2), (0, 4), (0, 6), (1, 3), (1, 5), (1, 7), (2, 4), (2, 6), (3, 5), (3, 7), (4, 6), (5, 7)]
        bondsun = []
        bonds_j2un = []


        self.energy_renorm = len(bonds) + len(bondsun) + len(bonds_j2) * j2 + len(bonds_j2un) * j2
        return ls.Operator(self.basis, ([ls.Interaction(operator * 2, bonds)] if len(bonds) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2 * 2, bonds_j2)] if len(bonds_j2) > 0 else []) + \
                                       ([ls.Interaction(operatorun * 2, bondsun)] if len(bondsun) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2un * 2, bonds_j2un)] if len(bonds_j2un) > 0 else [])), \
               ([[ls.Operator(self.basis, [ls.Interaction(operator, [bond])]), 2] for bond in bonds] if len(bonds) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2, [bond])]), j2 * 2.] for bond in bonds_j2] if len(bonds_j2) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operatorun, [bond])]), 2] for bond in bondsun] if len(bondsun) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2un, [bond])]), j2 * 2.] for bond in bonds_j2un] if len(bonds_j2un) > 0 else []), \
               bonds + bondsun + bonds_j2 + bonds_j2un



class HeisenbergHoneycomb_3x3(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = +1., j_zz = 1., j2=0., xBC='PBC', yBC='PBC'):
        assert ((Lx == 3) and (Ly == 3))
        operator = P_ij
        operator_j2 = P_ij
        n_sites = 18

        bonds = [(0, 1), (0, 13), (0, 15), (1, 6), (1, 10), (2, 3), (2, 15), (2, 17), (3, 6), (3, 8), (4, 5), (4, 13), (4, 17), (5, 8), (5, 10), (6, 7), (7, 12), (7, 16), (8, 9), (9, 12), (9, 14), (10, 11), (11, 14), (11, 16), (12, 13), (14, 15), (16, 17)]
        bonds_j2 = [(0, 2), (0, 6), (0, 14), (0, 12), (0, 4), (0, 10), (1, 3), (1, 15), (1, 13), (1, 5), (1, 7), (1, 11), (2, 14), (2, 16), (2, 4), (2, 8), (2, 6), (3, 15), (3, 17), (3, 5), (3, 9), (3, 7), (4, 8), (4, 10), (4, 12), (4, 16), (5, 17), (5, 13), (5, 11), (5, 9), (6, 8), (6, 12), (6, 16), (6, 10), (7, 9), (7, 13), (7, 17), (7, 11), (8, 10), (8, 14), (8, 12), (9, 11), (9, 15), (9, 13), (10, 14), (10, 16), (11, 15), (11, 17), (12, 14), (12, 16), (13, 15), (13, 17), (14, 16), (15, 17)]

        bondsun = []
        bonds_j2un = []


        self.energy_renorm = len(bonds) + len(bondsun) + len(bonds_j2) * j2 + len(bonds_j2un) * j2
        return ls.Operator(self.basis, ([ls.Interaction(operator * 2, bonds)] if len(bonds) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2 * 2, bonds_j2)] if len(bonds_j2) > 0 else []) + \
                                       ([ls.Interaction(operatorun * 2, bondsun)] if len(bondsun) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2un * 2, bonds_j2un)] if len(bonds_j2un) > 0 else [])), \
               ([[ls.Operator(self.basis, [ls.Interaction(operator, [bond])]), 2] for bond in bonds] if len(bonds) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2, [bond])]), j2 * 2.] for bond in bonds_j2] if len(bonds_j2) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operatorun, [bond])]), 2] for bond in bondsun] if len(bondsun) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2un, [bond])]), j2 * 2.] for bond in bonds_j2un] if len(bonds_j2un) > 0 else []), \
               bonds + bondsun + bonds_j2 + bonds_j2un, [2] * len(bonds) + [2 * j2] * len(bonds_j2)


class HeisenbergSquare_5x4(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = +1., j_zz = 1., j2=0., xBC='PBC', yBC='PBC'):
        #assert Lx % 2 == 0  # here we only ocnsider bipartite systems
        #assert Ly % 2 == 0
        assert Lx == 5 or Lx == 3
        assert Ly == 4

        operator = P_ij
        operator_j2 = P_ij
        operatorun = P_ijun
        operator_j2un = P_ijun

        n_sites = Lx * Ly

        bonds = []
        bonds_j2 = []
        bondsun = []
        bonds_j2un = []

        for site in range(n_sites):
            x, y = site % Lx, site // Lx

            site_up = ((x + 1) % Lx) + y * Lx
            site_right = x + ((y + 1) % Ly) * Lx

            if x + 1 < Lx:
                if self.unitary[site, site_up] == +1:
                    bonds.append((site, site_up))
                else:
                    bondsun.append((site, site_up))
            if y + 1 < Ly or yBC == 'PBC':
                if self.unitary[site, site_right] == +1:
                    bonds.append((site, site_right))
                else:
                    bondsun.append((site, site_right))


            if not np.isclose(j2, 0.0):
                site_up = ((x + 1) % Lx) + ((y + 1) % Ly) * Lx
                site_right = ((x + 1) % Lx) + ((y - 1) % Ly) * Lx
                if (x + 1 >= Lx):
                    continue

                if (y + 1 < Ly) or yBC == 'PBC':
                    if self.unitary[site, site_up] == +1:
                        bonds_j2.append((site, site_up))
                    else:
                        bonds_j2un.append((site, site_up))
                if (y - 1 >= 0) or yBC == 'PBC':
                    if self.unitary[site, site_right] == +1:
                        bonds_j2.append((site, site_right))
                    else:
                        bonds_j2un.append((site, site_right))

        print(bonds)
        print(bonds_j2)
        self.energy_renorm = len(bonds) + len(bondsun) + len(bonds_j2) * j2 + len(bonds_j2un) * j2
        return ls.Operator(self.basis, ([ls.Interaction(operator * 2, bonds)] if len(bonds) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2 * 2, bonds_j2)] if len(bonds_j2) > 0 else []) + \
                                       ([ls.Interaction(operatorun * 2, bondsun)] if len(bondsun) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2un * 2, bonds_j2un)] if len(bonds_j2un) > 0 else [])), \
               ([[ls.Operator(self.basis, [ls.Interaction(operator, [bond])]), 2] for bond in bonds] if len(bonds) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2, [bond])]), j2 * 2.] for bond in bonds_j2] if len(bonds_j2) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operatorun, [bond])]), 2] for bond in bondsun] if len(bondsun) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2un, [bond])]), j2 * 2.] for bond in bonds_j2un] if len(bonds_j2un) > 0 else []), \
               bonds + bondsun + bonds_j2 + bonds_j2un, \
               [2] * len(bonds) + [2 * j2] * len(bonds_j2)

class HeisenbergChain(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = +1., j_zz = 1., j2=0., xBC='PBC', yBC = 'PBC'):
        #assert Lx % 2 == 0  # here we only ocnsider bipartite systems
        #assert Ly % 2 == 0

        operator = P_ij
        operator_j2 = P_ij
        operatorun = P_ijun
        operator_j2un = P_ijun

        n_sites = Lx * Ly

        bonds = []
        bonds_j2 = []
        bondsun = []
        bonds_j2un = []

        for site in range(n_sites):
            x, y = site % Lx, site // Lx

            site_right = x + ((y + 1) % Ly) * Lx

            if y + 1 < Ly or yBC == 'PBC':
                if self.unitary[site, site_right] == +1:
                    bonds.append((site, site_right))
                else:
                    bondsun.append((site, site_right))


            if not np.isclose(j2, 0.0):
                site_up = ((x) % Lx) + ((y + 2) % Ly) * Lx
                if (y + 2 < Ly or yBC == 'PBC'):
                    if self.unitary[site, site_up] == +1:
                        bonds_j2.append((site, site_up))
                    else:
                        bonds_j2un.append((site, site_up))

        print(bonds + bonds_j2)
        self.energy_renorm = len(bonds) + len(bondsun) + len(bonds_j2) * j2 + len(bonds_j2un) * j2
        return ls.Operator(self.basis, ([ls.Interaction(operator * 2, bonds)] if len(bonds) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2 * 2, bonds_j2)] if len(bonds_j2) > 0 else []) + \
                                       ([ls.Interaction(operatorun * 2, bondsun)] if len(bondsun) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2un * 2, bonds_j2un)] if len(bonds_j2un) > 0 else [])), \
               ([[ls.Operator(self.basis, [ls.Interaction(operator, [bond])]), 2] for bond in bonds] if len(bonds) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2, [bond])]), j2 * 2.] for bond in bonds_j2] if len(bonds_j2) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operatorun, [bond])]), 2] for bond in bondsun] if len(bondsun) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2un, [bond])]), j2 * 2.] for bond in bonds_j2un] if len(bonds_j2un) > 0 else []), \
               bonds + bondsun + bonds_j2 + bonds_j2un, \
               [2] * len(bonds) + [2 * j2] * len(bonds_j2)


class TFIMChain(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, h, xBC='PBC', yBC = 'PBC'):
        operator = np.kron(sz, sz)
        operator_h = sx

        n_sites = Lx * Ly

        bonds = []
        bonds_j2 = []


        bonds_j2 = list(np.arange(Lx * Ly))
        for site in range(n_sites):
            bonds.append((site, (site + 1) % n_sites))

        self.energy_renorm = 0.0
        return ls.Operator(self.basis, ([ls.Interaction(operator, bonds)] if len(bonds) > 0 else []) + \
                                       ([ls.Interaction(h * operator_h, bonds_j2)] if len(bonds_j2) > 0 else [])), \
               ([[ls.Operator(self.basis, [ls.Interaction(operator, [bond])]), 1] for bond in bonds] if len(bonds) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_h, [bond])]), h] for bond in bonds_j2] if len(bonds_j2) > 0 else []) + [] + [], \
               bonds + bonds_j2, \
               [1] * len(bonds) + [h] * len(bonds_j2)


class HeisenbergKagome18(Hamiltonian):
    def _get_Hamiltonian_matrix(self, Lx, Ly, j_pm = +1., j_zz = 1., j2=0., xBC='PBC', yBC = 'PBC'):
        #assert Lx % 2 == 0  # here we only ocnsider bipartite systems
        #assert Ly % 2 == 0

        operator = P_ij
        operator_j2 = P_ij
        operatorun = P_ijun
        operator_j2un = P_ijun

        n_sites = Lx * Ly

        bonds = [(0, 6), (6, 1), (1, 7), (7, 2), (2, 8), (8, 0), (3, 9), (9, 4), (4, 10), (10, 5), (5, 11), (11, 3), (0, 12), (12, 3), (3, 15), (15, 0), (1, 13), (13, 4), (4, 16), (16, 1), (2, 14), (14, 5), (5, 17), (17, 2), (6, 12), (12, 11), (11, 17), (17, 7), (7, 13), (13, 9), (9, 15), (15, 8), (8, 14), (14, 10), (10, 16), (16, 6)]

        bonds_j2 = [(0, 1), (6, 7), (1, 2), (7, 8), (2, 0), (8, 6), (3, 4), (9, 10), (4, 5), (10, 11), (5, 3), (11, 9), (0, 3), (12, 15), (1, 4), (13, 16), (2, 5), (14, 17), (6, 11), (12, 17), (11, 7), (17, 13), (7, 9), (13, 15), (9, 8), (14, 15), (8, 10), (14, 16), (10, 6), (16, 12)]


        bondsun = []
        bonds_j2un = []

        self.energy_renorm = len(bonds) + len(bondsun) + len(bonds_j2) * j2 + len(bonds_j2un) * j2
        return ls.Operator(self.basis, ([ls.Interaction(operator * 2, bonds)] if len(bonds) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2 * 2, bonds_j2)] if len(bonds_j2) > 0 else []) + \
                                       ([ls.Interaction(operatorun * 2, bondsun)] if len(bondsun) > 0 else []) + \
                                       ([ls.Interaction(j2 * operator_j2un * 2, bonds_j2un)] if len(bonds_j2un) > 0 else [])), \
               ([[ls.Operator(self.basis, [ls.Interaction(operator, [bond])]), 2] for bond in bonds] if len(bonds) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2, [bond])]), j2 * 2.] for bond in bonds_j2] if len(bonds_j2) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operatorun, [bond])]), 2] for bond in bondsun] if len(bondsun) > 0 else []) + \
               ([[ls.Operator(self.basis, [ls.Interaction(operator_j2un, [bond])]), j2 * 2.] for bond in bonds_j2un] if len(bonds_j2un) > 0 else []), \
               bonds + bondsun + bonds_j2 + bonds_j2un, \
               [2] * len(bonds) + [2 * j2] * len(bonds_j2)

