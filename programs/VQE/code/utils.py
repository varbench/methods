import os
import sys
import numpy as np
from time import time
import lattice_symmetries as ls
#import qiskit
import scipy
import scipy.linalg
from numba import jit

@jit(nopython=True)
def get_trans(idxs_mapped, idxs_mapped_flip, st_list):
    permutation = [0]
    for idx, idx_flip in zip(idxs_mapped, idxs_mapped_flip):
        if idx in st_list:
            permutation.append(idx)
        else:
            permutation.append(idx_flip)
    return permutation[1:]


def index_to_spin(index, number_spins = 16):
    return (((index.reshape(-1, 1).astype(np.int64) & (1 << np.arange(number_spins).astype(np.int64)))) > 0)

def spin_to_index(spin, number_spins = 16):
    a = 2 ** np.arange(number_spins)
    return spin.dot(a)

def last_to_ij(i, j, number_spins = 16):
    idxs = np.arange(2 ** number_spins)
    spin = index_to_spin(idxs, number_spins)
    spin[:, i], spin[:, N - 2] = spin[:, N - 2], spin[:, i]
    spin[:, j], spin[:, N - 1] = spin[:, N - 1], spin[:, j]
    return spin_to_index(spin, number_spins)

def import_config(filename: str):
    import importlib

    module_name, extension = os.path.splitext(os.path.basename(filename))
    module_dir = os.path.dirname(filename)
    if extension != ".py":
        raise ValueError(
            "Could not import the module from {!r}: not a Python source file.".format(
                filename
            )
        )
    if not os.path.exists(filename):
        raise ValueError(
            "Could not import the module from {!r}: no such file or directory".format(
                filename
            )
        )
    sys.path.insert(0, module_dir)
    module = importlib.import_module(module_name)
    sys.path.pop(0)
    return module


def get_x_symmetry_map(Lx, Ly, basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = Lx * Ly)

    map_site = []
    for i in range(Lx * Ly):
        x, y = i % Lx, i // Lx
        j = (x + 1) % Lx + y * Lx
        map_site.append(j)
    map_site = np.array(map_site)
    print('x', map_site)

    spins = spins[:, map_site]

    idxs_mapped = spin_to_index(spins, number_spins = Lx * Ly).astype(np.uint64)
    idxs_mapped2 = np.zeros((len(idxs_mapped), 8), dtype=np.uint64)
    idxs_mapped2[:, 0] = idxs_mapped
    rep, _, _ = basis.batched_state_info(idxs_mapped2)
    rep = rep[:, 0]
    assert len(np.unique(rep)) == len(rep)
    return map_site, np.argsort(rep)
    #idxs_mapped_flip = spin_to_index(1 - spins, number_spins = Lx * Ly).astype(np.uint64)

    #permutation = [] #get_trans(list(idxs_mapped), list(idxs_mapped_flip), list(basis.states.astype(np.uint64)))
    #st_list = basis.states
    #for idx, idx_flip in zip(idxs_mapped, idxs_mapped_flip):
    #    if idx in st_list:
    #        permutation.append(idx)
    #    else:
    #        permutation.append(idx_flip)

    return map_site, np.argsort(np.array(permutation)) #np.argsort(spin_to_index(spins, number_spins = Lx * Ly))

def get_y_symmetry_map(Lx, Ly, basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = Lx * Ly)

    map_site = []
    for i in range(Lx * Ly):
        x, y = i % Lx, i // Lx
        j = x + ((y + 1) % Ly) * Lx
        map_site.append(j)
    map_site = np.array(map_site)
    print('y', map_site)

    spins = spins[:, map_site]

    idxs_mapped = spin_to_index(spins, number_spins = Lx * Ly).astype(np.uint64)
    idxs_mapped2 = np.zeros((len(idxs_mapped), 8), dtype=np.uint64)
    idxs_mapped2[:, 0] = idxs_mapped
    rep, _, _ = basis.batched_state_info(idxs_mapped2)
    rep = rep[:, 0]
    assert len(np.unique(rep)) == len(rep)
    return map_site, np.argsort(rep)

    idxs_mapped = spin_to_index(spins, number_spins = Lx * Ly)
    idxs_mapped_flip = spin_to_index(1 - spins, number_spins = Lx * Ly)

    permutation = []
    st_list = basis.states
    for idx, idx_flip in zip(idxs_mapped, idxs_mapped_flip):
        if idx in st_list:
            permutation.append(idx)
        else:
            permutation.append(idx_flip)

    return map_site, np.argsort(np.array(permutation))


def get_Cx_symmetry_map(Lx, Ly, basis, su2=False):
    n_qubits = Lx * Ly
    spins = index_to_spin(basis.states, number_spins = Lx * Ly)

    map_site = []
    for i in range(Lx * Ly):
        x, y = i % Lx, i // Lx
        j = (Lx - 1 - x) % Lx + y * Lx
        map_site.append(j)
    map_site = np.array(map_site)
    print('Cx', map_site)

    spins = spins[:, map_site]

    idxs_mapped = spin_to_index(spins, number_spins = Lx * Ly).astype(np.uint64)
    idxs_mapped2 = np.zeros((len(idxs_mapped), 8), dtype=np.uint64)
    idxs_mapped2[:, 0] = idxs_mapped
    rep, _, _ = basis.batched_state_info(idxs_mapped2)
    rep = rep[:, 0]
    assert len(np.unique(rep)) == len(rep)
    return map_site, np.argsort(rep)

    idxs_mapped = spin_to_index(spins, number_spins = Lx * Ly)
    idxs_mapped_flip = spin_to_index(1 - spins, number_spins = Lx * Ly)

    permutation = []
    st_list = basis.states
    for idx, idx_flip in zip(idxs_mapped, idxs_mapped_flip):
        if idx in st_list:
            permutation.append(idx)
        else:
            permutation.append(idx_flip)

    return map_site, np.argsort(np.array(permutation))

def r_to_index(x, y, Lx, Ly):
    return y * Lx + x

def index_to_r(idx, Lx, Ly):
    return idx % Lx, idx // Lx


def get_rot_symmetry_map(Lx, Ly, basis, su2=False):
    assert Lx == Ly
    assert Lx == 4
    '''
    xmap = []
    L = Lx
    for idx in range(L ** 2):
        x, y = index_to_r(idx, L, L)
        xpr = (-y) % L
        ypr = x

        xmap.append(r_to_index(xpr, ypr, L, L))

    
    xmap = np.array(xmap)
    '''

    xmap = np.array([12, 8, 4, 0, 13, 9, 5, 1, 14, 10, 6, 2, 15, 11, 7, 3])
    n_qubits = Lx * Ly
    spins = index_to_spin(basis.states, number_spins = Lx * Ly)

    spins = spins[:, xmap]

    idxs_mapped = spin_to_index(spins, number_spins = Lx * Ly).astype(np.uint64)
    idxs_mapped2 = np.zeros((len(idxs_mapped), 8), dtype=np.uint64)
    idxs_mapped2[:, 0] = idxs_mapped
    rep, _, _ = basis.batched_state_info(idxs_mapped2)
    rep = rep[:, 0]
    assert len(np.unique(rep)) == len(rep)
    return xmap, np.argsort(rep)


def get_Cy_symmetry_map(Lx, Ly, basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = Lx * Ly)

    map_site = []
    for i in range(Lx * Ly):
        x, y = i % Lx, i // Lx
        j = x % Lx + ((Ly - 1 - y) % Ly) * Lx
        map_site.append(j)
    map_site = np.array(map_site)

    spins = spins[:, map_site]

    idxs_mapped = spin_to_index(spins, number_spins = Lx * Ly).astype(np.uint64)
    idxs_mapped2 = np.zeros((len(idxs_mapped), 8), dtype=np.uint64)
    idxs_mapped2[:, 0] = idxs_mapped
    rep, _, _ = basis.batched_state_info(idxs_mapped2)
    rep = rep[:, 0]
    assert len(np.unique(rep)) == len(rep)
    print('Cy', map_site)
    return map_site, np.argsort(rep)


    idxs_mapped = spin_to_index(spins, number_spins = Lx * Ly)
    idxs_mapped_flip = spin_to_index(1 - spins, number_spins = Lx * Ly)

    permutation = []
    st_list = basis.states
    for idx, idx_flip in zip(idxs_mapped, idxs_mapped_flip):
        if idx in st_list:
            permutation.append(idx)
        else:
            permutation.append(idx_flip)

    print('Cy', map_site)
    return map_site, np.argsort(np.array(permutation))

def get_hexagon_rot_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 6)

    map_site = np.array([1, 2, 3, 4, 5, 0])
    assert np.allclose(map_site[map_site][map_site][map_site][map_site][map_site], np.arange(6))
  

    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 6))

def get_hexagon_mir_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 6)

    map_site = np.array([0, 5, 4, 3, 2, 1])
    assert np.allclose(map_site[map_site], np.arange(6))


    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 6))


def get_honeycomb_2x2_trx_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 8)

    map_site = np.array([4, 5, 6, 7, 0, 1, 2, 3])
    assert np.allclose(map_site[map_site], np.arange(8))


    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 8))


def get_honeycomb_3x3_trx_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 18)

    map_site = np.array([2, 3, 4, 5, 0, 1, 8, 9, 10, 11, 6, 7, 14, 15, 16, 17, 12, 13])
    assert np.allclose(map_site[map_site][map_site], np.arange(18))


    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 18))


def get_honeycomb_2x2_try_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 8)

    map_site = np.array([2, 3, 0, 1, 6, 7, 4, 5])
    assert np.allclose(map_site[map_site], np.arange(8))


    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 8))


def get_honeycomb_3x3_try_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 18)

    map_site = np.array([6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5])
    assert np.allclose(map_site[map_site][map_site], np.arange(18))


    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 18))



def get_honeycomb_2x2_rot_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 8)

    map_site = np.array([1, 4, 7, 2, 3, 6, 5, 0])
    assert np.allclose(map_site[map_site][map_site][map_site][map_site][map_site], np.arange(8))


    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 8))


def get_honeycomb_3x3_rot_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 18)

    map_site = np.array([1, 6, 15, 2, 11, 16, 3, 8, 17, 4, 7, 12, 5, 10, 13, 0, 9, 14])
    assert np.allclose(map_site[map_site][map_site][map_site][map_site][map_site], np.arange(18))


    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 18))


def get_honeycomb_2x2_mir1_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 8)

    map_site = np.array([0, 5, 6, 3, 4, 1, 2, 7])
    assert np.allclose(map_site[map_site], np.arange(8))


    spins = spins[:, map_site] 
    return map_site, np.argsort(spin_to_index(spins, number_spins = 8))


def get_honeycomb_3x3_mir2_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 18)

    map_site = np.array([0, 13, 14, 9, 10, 5, 12, 7, 8, 3, 4, 17, 6, 1, 2, 15, 16, 11])
    assert np.allclose(map_site[map_site], np.arange(18))


    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 18))

def get_honeycomb_3x3_mir1_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 18)

    map_site = np.array([17, 16, 13, 12, 15, 14, 7, 6, 9, 8, 11, 10, 3, 2, 5, 4, 1, 0])
    assert np.allclose(map_site[map_site], np.arange(18))


    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 18))


def get_honeycomb_2x2_mir2_symmetry_map(basis, su2=False):
    spins = index_to_spin(basis.states, number_spins = 8)

    map_site = np.array([7, 2, 1, 4, 3, 6, 5, 0])
    assert np.allclose(map_site[map_site], np.arange(8))


    spins = spins[:, map_site]
    return map_site, np.argsort(spin_to_index(spins, number_spins = 8))


def compute_norm_sample(state, projector, N_samples, noise_p = 0.):
    t = time()
    ps = [(0.5 - np.vdot(projector(state, proj_idx), state) / 2.).real for proj_idx in range(projector.nterms)]
    ps = np.clip(ps, a_min=0., a_max=1.)
    N_ups = np.random.binomial(N_samples, ps)
    #print('sample estimation of N(theta) = ', time() - t)

    flip_mask = np.ones(projector.nterms) if noise_p == 0. else 1. - 2 * (np.random.uniform(0, 1, size=projector.nterms) < noise_p)

    return np.sum(flip_mask * (1 - 2 * N_ups / N_samples)) / projector.nterms

def compute_overlap_sample(state1, state2, N_samples):
    t = time()
    p = np.clip(np.abs(np.vdot(state1, state2)) ** 2, a_min=0., a_max=1.)
    return 1. * np.random.binomial(N_samples, p) / N_samples

'''
def compute_norm_sample_noisy(states, projector, N_samples, noise_p, basis_with_ancilla, noise_operators):
    t = time()
    norms = []
    states_ancilla = np.zeros(states.shape[0], 2 * states.shape[1], dtype=np.complex128)

    states_ancilla[:, :len(states.shape[1])] = states / np.sqrt(2)
    states_ancilla[:, len(states.shape[1]):] = states / np.sqrt(2)

    proj_ctr = 0

    for actions, vector in zip(projector.all_double_permutations, states_ancilla):
        N = len(actions)
        site = np.random.choose(np.arange(N), p = [(1 - noise_p) ** N] + [(1 - (1 - noise_p) ** N) / (N - 1)] * (N - 1))
        if site == 0:
            states_ancilla[proj_ctr, len(states.shape[1]):] = states_ancilla[proj_ctr, len(states.shape[1]):][actions[0][1]]  # act after
            continue

        ### then we act noise ###
        states_ancilla[proj_ctr, len(states.shape[1]):] = states_ancilla[proj_ctr, len(states.shape[1]):][actions[0][0]]  # act before
        pair = projector.all_pair_permutations[proj_ctr][site]

        sites = tuple(np.random.choose([0, pair[0], pair[1], 2, p=[1. / 3., 1. / 3., 1. / 3.]]))
        op = ls.Operator(basis_with_ancilla, [ls.Interaction(noise_operators[np.random.choose(len(states_ancilla))], [sites])]) # get noise



        states_ancilla = op(states_ancilla)  # act noise
        states_ancilla[proj_ctr, len(states.shape[1]):] = states_ancilla[proj_ctr, len(states.shape[1]):][actions[0][1]]  # after
        proj_ctr += 1



    ps = [(0.5 - np.vdot(projector(state, proj_idx), state) / 2.).real for proj_idx in range(projector.nterms)]
    ps = np.clip(ps, a_min=0., a_max=1.)
    N_ups = np.random.binomial(N_samples, ps)
    print('sample estimation of N(theta) = ', time() - t)
    return np.sum(1 - 2 * N_ups / N_samples) / projector.nterms

    for proj_idx in range(projector.nterms):
        state_proj = projector(state, proj_idx)
        state_ancilla[:len(state)] = (state + state_proj) / 2.
        state_ancilla[len(state):] = (state - state_proj) / 2.


        idxs = np.random.choice(indexes, p=np.abs(state_ancilla) ** 2, replace=True, size=N_samples)

        N_up = np.sum(idxs > len(state))

        norms.append(1 - 2 * N_up / N_samples)

    print('sample estimation of N(theta) = ', time() - t)
    return np.sum(norms) / projector.nterms
'''


from copy import deepcopy
def get_symmetry_unique_bonds(bonds, permutations):
    removed_bonds = []
    groups = []
    idxs_groups = []
    tot_size = 0

    for bond in bonds:
        group = []
        idxs_group = []
        bonds_symm = [(perm[bond[0]], perm[bond[1]]) for perm in permutations]
        
        for bond_symm in bonds_symm:
            if bond_symm in removed_bonds or bond_symm[::-1] in removed_bonds:
                continue

            if bond_symm in bonds:
                removed_bonds.append(bond_symm)
                group.append(bond_symm)
                idxs_group.append(bonds.index(bond_symm))
            else:
                removed_bonds.append(bond_symm[::-1])
                group.append(bond_symm[::-1])
                idxs_group.append(bonds.index(bond_symm[::-1]))
        if len(group) == 0:
            continue
        groups.append(deepcopy(group))
        idxs_groups.append(deepcopy(idxs_group))
        tot_size += len(group)
    assert tot_size == len(bonds)
    return groups, idxs_groups




def compute_energy_sample(state, hamiltonian, projector, N_samples, noise_p = 0.):
    t = time()

    state0_proj_inv = np.array([projector(state, proj_idx, inv=True) for proj_idx in range(projector.nterms)]).conj()

    js = np.array([hamiltonian(None, ham_idx, vector=False) for ham_idx in range(hamiltonian.nterms)])
    states_ham = np.array([hamiltonian(state, ham_idx) for ham_idx in range(hamiltonian.nterms)])

    ps = (0.5 - np.dot(state0_proj_inv, states_ham.T) / 2.).real
    flip_mask = np.ones(shape=(projector.nterms, len(js))) if noise_p == 0. else 1. - 2 * (np.random.uniform(0, 1, size=(projector.nterms, len(js))) < noise_p)

    ps = np.clip(ps, a_min=0., a_max=1.)
    N_ups = np.random.binomial(N_samples, ps)
    energy = np.sum((flip_mask * (1 - 2 * N_ups / N_samples)).dot(js))
    #print('sample estimation of E(theta) = ', time() - t)
    return energy / projector.nterms


def get_hamiltonian_power_expectation_exact(H_power, circuit, hamiltonian, projector, config):
    state = circuit.__call__()
    state_proj = projector(state)

    for _ in range(H_power):
        state_proj = hamiltonian(state_proj)# - hamiltonian.energy_renorm * state_proj
    return np.vdot(state, state_proj)


def get_hamiltonian_power_expectation_sampling(H_power, circuit, hamiltonian, projector, config):
    if H_power == 0:
        return compute_norm_sample(circuit.__call__(), projector, config.N_samples)

    state = circuit.__call__()
    state0_proj_inv = np.array([projector(state, proj_idx, inv=True) for proj_idx in range(projector.nterms)]).conj()

    import itertools    

    n_stored_powers = 3
    stored_wfs = []
    stored_js = []

    H_power_total = 0.0
    max_idx =  np.min([n_stored_powers, H_power])

    precomputed_wfs = [state]
    for _ in range(max_idx):
        precomputed_wfs.append(hamiltonian(precomputed_wfs[-1], 0))
    idxs = np.zeros(max_idx, dtype=np.int64)


    for perm_combination in itertools.product(*([np.arange(len(hamiltonian.bonds))] * max_idx)):
        #print(perm_combination)

        if len(np.where(idxs != perm_combination)[0]) != 0:
            recompute_from = np.min(np.where(idxs != perm_combination)[0])
            idxs = perm_combination

            for recompute_idx in range(recompute_from, max_idx):
                precomputed_wfs[1 + recompute_idx] = hamiltonian(precomputed_wfs[recompute_idx], perm_combination[recompute_idx])
        state_act = precomputed_wfs[-1]


        #state_act = state.copy()
        j_total = 1.0
        for perm in perm_combination:
            #state_act = hamiltonian(state_act, perm)
            j_total *= hamiltonian(None, perm, vector=False)

        ps = 0.5 - np.dot(state0_proj_inv, state_act).real / 2.
        N_ups = np.random.binomial(config.N_samples, ps)


        if H_power <= n_stored_powers:
            H_power_total += np.sum((1 - 2 * N_ups / config.N_samples)) * j_total / projector.nterms
        else:
            stored_wfs.append(state_act.copy())
            stored_js.append(j_total)

    if H_power <= n_stored_powers:
        return H_power_total

    for perm_combination in itertools.product(*([np.arange(len(hamiltonian.bonds))] * (H_power - n_stored_powers))):
        #print(perm_combination)
        states_act = np.array(stored_wfs).copy() * 1.0
        js_act = np.array(stored_js).copy() * 1.0

        for perm in perm_combination:
            states_act = hamiltonian(states_act.T, perm).T
            js_act *= hamiltonian(None, perm, vector=False)
        ps = 0.5 - np.dot(state0_proj_inv, states_act.T).real / 2.
        N_ups = np.random.binomial(config.N_samples, ps)

        H_power_total += np.sum((1 - 2 * N_ups / config.N_samples) * js_act) / projector.nterms
    return H_power_total



def compute_energy_sample_symmetrized(state, hamiltonian, projector, N_samples):
    bond_groups, idxs_groups = get_symmetry_unique_bonds(hamiltonian.bonds, projector.maps)

    t = time()
    energies = []
    state_ancilla = np.zeros(2 * len(state), dtype=np.complex128)
    indexes = np.arange(len(state_ancilla))
    
    for idxs_group in idxs_groups:
        ham_idx = idxs_group[0]
        state_ham, j = hamiltonian(state, ham_idx)
        for proj_idx in range(projector.nterms):
            state_proj = projector(state_ham, proj_idx)
            state_ancilla[:len(state)] = (state + state_proj) / 2.
            state_ancilla[len(state):] = (state - state_proj) / 2.

            idxs = np.random.choice(indexes, p=np.abs(state_ancilla) ** 2, replace=True, size=N_samples * 5)
  
            N_up = np.sum(idxs >= len(state))

            energies.append((1 - 2 * N_up / N_samples / 5.) * j * len(idxs_group))

    print('sample estimation of E(theta) = ', time() - t)
    return np.sum(energies) / projector.nterms


def compute_metric_tensor_sample(states, projector, N_samples, noise_p = 0.):
    t = time()
    t_proj = 0.
    t_norm = 0.
    t_samp = 0.

    ctr = 0
    ps_real = np.empty((len(states), len(states), projector.nterms), dtype=np.float64)
    ps_imag = np.empty((len(states), len(states), projector.nterms), dtype=np.float64)
    #states = np.asfortranarray(states)

    statesL = np.zeros((len(projector.lpermutations), states.shape[0], states.shape[1]), dtype=np.complex64)
    statesR = np.ascontiguousarray(np.zeros((len(projector.rpermutations), states.shape[0], states.shape[1]), dtype=np.complex64))

    states = np.ascontiguousarray(states)

    for idxl, perm in enumerate(projector.lpermutations):
        z = time()
        for j in range(states.shape[0]):
            statesL[idxl, j, :] = states[j, perm].astype(np.complex64)
        t_proj += time() - z
    z = time()
    for idxr, perm in enumerate(projector.rpermutations):
        for j in range(states.shape[0]):
            statesR[idxr, j, :] = states[j, perm].astype(np.complex64).conj()
    t_proj += time() - z

    #for proj_idx in range(projector.nterms):
    for pair_idxs in projector.list_of_pairs:
        #z = time()
        #states_j_proj = projector(states, proj_idx)#.conj() # np.array([projector(states[j], proj_idx) for j in range(len(states))]).conj(
        #t_proj += time() - z

        z = time()
        proj_idx = pair_idxs[0]
        l, r = projector.left_right_decompositions[proj_idx]
        
        #ps = 0.5 - scipy.linalg.blas.zher2k(1., states, states_j_proj).conj() / 4.
        ps = 0.5 - np.dot(statesL[l], statesR[r].T).conj() / 2.
        #ps = 0.5 - np.dot(states, states_j_proj.T.conj()).conj() / 2.
        #ps = 0.5 - cp.asnumpy(cp.asarray(states).dot(cp.asarray(states_j_proj.conj().T))).conj() / 2.
        #ps = 0.5 - scipy.linalg.blas.zgemm(1, np.asfortranarray(states), np.asfortranarray(states_j_proj), trans_b = 2).conj() / 2.
        ps_real[..., proj_idx] = ps.real
        ps_imag[..., proj_idx] = ps.imag + 0.5

        if len(pair_idxs) == 2:
            ps_real[..., pair_idxs[1]] = ps.real.T
            ps_imag[..., pair_idxs[1]] = -(ps.imag + 0.5).T

        t_norm += time() - z

    z = time()
    N_ups_reals = np.random.binomial(N_samples, np.clip(ps_real, a_min=0., a_max=1.))
    N_ups_imags = np.random.binomial(N_samples, np.clip(ps_imag, a_min=0., a_max=1.))
    t_samp += time() - z

    flip_mask = np.ones(shape=(states.shape[0], states.shape[0], projector.nterms)) if noise_p == 0. else 1. - 2 * (np.random.uniform(0, 1, size=(states.shape[0], states.shape[0], projector.nterms)) < noise_p)

    #print('flip mask mean for connectivity: ', np.mean(flip_mask), noise_p)
    MT = (flip_mask * (1 - 2 * N_ups_reals / N_samples)).mean(axis=-1) + 1.0j * (flip_mask * (1 - 2 * N_ups_imags / N_samples)).mean(axis=-1)
    #print('sample estimation of MT(theta) = ', time() - t)
    #print(t_proj, t_norm, t_samp)

    #print(np.linalg.norm(MT.real - MT.real.conj().T))
    #exit(-1)
    return MT


def compute_connectivity_sample(state0, states, projector, N_samples, noise_p = 0.):
    t = time()
    # thetafull = np.exp(1.0j * np.pi / 2. * theta)

    state0_proj_inv = np.array([projector(state0, proj_idx, inv=True) for proj_idx in range(projector.nterms)]).conj()


    flip_mask = np.ones(shape=(projector.nterms, states.shape[0])) if noise_p == 0. else 1. - 2 * (np.random.uniform(0, 1, size=(projector.nterms, states.shape[0])) < noise_p)


    #print('flip mask mean for connectivity: ', np.mean(flip_mask))
    ps_general = np.dot(state0_proj_inv, states.T)
    ps_real = (0.5 - ps_general / 2.).real
    ps_imag = (0.5 + 1.0j * ps_general / 2.).real

    ps_real = np.clip(ps_real, a_min=0., a_max=1.)
    ps_imag = np.clip(ps_imag, a_min=0., a_max=1.)


    #ps = (0.5 - np.dot(state0_proj_inv, states.T) / 2.).real
    #ps = np.clip(ps, a_min=0., a_max=1.)
    N_ups_real = np.random.binomial(N_samples, ps_real)
    N_ups_imag = np.random.binomial(N_samples, ps_imag)
    connectivity = np.sum(flip_mask * (1 - 2 * N_ups_real / N_samples), axis = 0) + 1.0j * np.sum(flip_mask * (1 - 2 * N_ups_imag / N_samples), axis = 0)

    #print('sample estimation of connectivity(theta) = ', time() - t)

    return connectivity / projector.nterms

def compute_energy_der_sample(state0, states, hamiltonian, projector, N_samples, noise_p = 0.):
    t = time()
    time_sampling = 0.
    t_proj = 0.
    t_ham = 0.
    t_samp = 0.
    t_norm = 0.

    z = time()
    state0_proj_inv = np.array([projector(state0, proj_idx, inv=True) for proj_idx in range(projector.nterms)]).conj()
    t_proj += time() - z


    ps = np.empty((projector.nterms, states.shape[1], hamiltonian.nterms), dtype=np.float64)
    js = []
    for ham_idx in range(hamiltonian.nterms):
        z = time()
        states_ham = hamiltonian(states, ham_idx) # np.empty((hamiltonian.nterms, len(state0)), dtype=np.complex128)
        j = hamiltonian(None, ham_idx, vector=False)
        js.append(j)
        t_ham += time() - z

        z = time()
        ps[..., ham_idx] = 0.5 - np.dot(state0_proj_inv, states_ham).real / 2.
        t_norm += time() - z

    z = time()
    N_ups = np.random.binomial(N_samples, ps)
    t_samp += time() - z

    flip_mask = np.ones(shape=(projector.nterms, states.shape[1], len(js))) if noise_p == 0. else 1. - 2 * (np.random.uniform(0, 1, size=(projector.nterms, states.shape[1], len(js))) < noise_p)

    #print('flip mask mean for der energy: ', np.mean(flip_mask), noise_p)
    der = np.sum((flip_mask * (1 - 2 * N_ups / N_samples)).dot(js), axis=0)

    #print(t_proj, t_ham, t_samp, t_norm)
    #print('sample estimation of energy der numerator(theta) = ', time() - t, 'sampling:', time_sampling)
    return der / projector.nterms


def compute_norm_qiskit_survivalrate(circuit, projector, N_samples, noise_model):
    basis_gates = noise_model.basis_gates

    circ = circuit.init_circuit_qiskit()
    circ = circuit.act_dimerization_qiskit(circ)
    circ = circuit.act_psi_qiskit(circ, circuit.params)

    projector_circuits = []
    for idx, cycl in enumerate(projector.cycl):
        current_circ = circ.copy('circuit_projector_{:d}'.format(idx))

        pair_permutations = projector.all_pair_permutations[idx]
        current_circ = circuit.act_permutation_qiskit(current_circ, pair_permutations)

        current_circ = circuit.act_psi_qiskit(current_circ, circuit.params, inverse=True)
        current_circ = circuit.act_dimerization_qiskit(current_circ, inverse=True)
        current_circ.save_amplitudes([0])

        projector_circuits.append(current_circ.copy())

    t = time()

    backend = qiskit.providers.aer.QasmSimulator(method='statevector', noise_model=noise_model)
    result = qiskit.execute(projector_circuits, backend, shots=1, basis_gates = backend.configuration().basis_gates).result()
    print('noisy simulation with quasm took', time() - t)

    survivals = []
    for i in range(len(projector.cycl)):
        survivals.append(np.abs(result.data(i)['amplitudes'][0]) ** 2)

    survivals = np.clip(survivals, a_min=0., a_max=1.)
    survivals = np.random.binomial(N_samples, np.sqrt(np.array(survivals))) / N_samples
    

    return np.mean(survivals)


def compute_norm_qiskit_hadamardtest(circuit, projector, N_samples, noise_model):
    circ = circuit.init_circuit_qiskit(ancilla=True)
    circ.h(circuit.n_qubits)
    circ = circuit.act_dimerization_qiskit(circ)
    circ = circuit.act_psi_qiskit(circ, circuit.params)

    projector_circuits = []
    for idx, cycl in enumerate(projector.cycl):
        current_circ = circ.copy('circuit_projector_{:d}'.format(idx))

        pair_permutations = projector.all_pair_permutations[idx]
        current_circ = circuit.act_permutation_qiskit(current_circ, pair_permutations, ancilla=True, ancilla_qubit_idx = circuit.n_qubits)
        current_circ.h(circuit.n_qubits)
        current_circ.measure(circuit.n_qubits, 0)
        projector_circuits.append(current_circ.copy())

    t = time()

    backend = qiskit.providers.aer.QasmSimulator(method='statevector', noise_model=noise_model, max_parallel_experiments=len(projector_circuits))
    result = qiskit.execute(projector_circuits, backend, shots=N_samples, basis_gates = backend.configuration().basis_gates).result()

    print('noisy simulation of norm with quasm took', time() - t)

    survivals = []
    for i in range(len(projector_circuits)):
        survivals.append(-1. + 2. * result.data(i)['counts']['0x0'] / N_samples if '0x0' in result.data(i)['counts'].keys() else -1)
    return np.mean(survivals * np.array(projector.characters))



def compute_connectivity_qiskit_hadamardtest(circuit, projector, N_samples, noise_model, real_imag, idx_from, idx_to):
    backend = qiskit.providers.aer.QasmSimulator(method='statevector', noise_model=noise_model, max_parallel_experiments=len(projector.cycl) * len(circuit.params))

    connectivities = []
    phase = 0 if real_imag == 'real' else -np.pi / 2.
    projector_circuits = []
    for i_param in range(len(circuit.params)):
        #print(i_param)
        circ = circuit.init_circuit_qiskit(ancilla=True)
        circ.h(circuit.n_qubits)
        circ = circuit.act_dimerization_qiskit(circ)
        circ = circuit.act_psi_qiskit_string(circ, circuit.params, 0, i_param)
        circ = circuit.act_derivative_qiskit(circ, i_param, ancilla_qubit_idx = circuit.n_qubits)
        circ = circuit.act_psi_qiskit_string(circ, circuit.params, i_param, len(circuit.params))
        circ.u1(phase, circuit.n_qubits)  # for real/imag

        for idx, cycl in enumerate(projector.cycl):
            if idx not in range(idx_from, idx_to):
                continue
            current_circ = circ.copy('circuit_projector_{:d}'.format(idx))

            pair_permutations = projector.all_pair_permutations[idx]
            current_circ = circuit.act_permutation_qiskit(current_circ, pair_permutations, ancilla=True, ancilla_qubit_idx = circuit.n_qubits)
            current_circ.h(circuit.n_qubits)

            current_circ.measure(circuit.n_qubits, 0)
            projector_circuits.append(current_circ.copy())

    t = time()
    result = qiskit.execute(projector_circuits, backend, shots=N_samples).result()
    print('noisy simulation of all connectivities with quasm took', time() - t)

    survivals = []

    for i in range(len(projector_circuits)):
        survivals.append(-1. + 2. * result.data(i)['counts']['0x0'] / N_samples if '0x0' in result.data(i)['counts'].keys() else -1)
    connectivities = np.array(survivals).reshape((len(circuit.params), -1))
    connectivities = np.einsum('ij,j->ij', connectivities, np.array(projector.characters)).mean(axis=-1)

    return -connectivities


def compute_der_qiskit_hadamardtest(circuit, hamiltonian, projector, N_samples, noise_model, idx_from, idx_to):
    backend = qiskit.providers.aer.QasmSimulator(method='statevector', noise_model=noise_model, \
                                                 max_parallel_experiments=len(projector.cycl) * len(hamiltonian.bonds))# * len(circuit.params))

    der = []
    for i_param in range(len(circuit.params)):
        #print(i_param)
        circ = circuit.init_circuit_qiskit(ancilla=True)
        circ.h(circuit.n_qubits)
        circ = circuit.act_dimerization_qiskit(circ)
        circ = circuit.act_psi_qiskit_string(circ, circuit.params, 0, i_param)
        circ = circuit.act_derivative_qiskit(circ, i_param, ancilla_qubit_idx = circuit.n_qubits)
        circ = circuit.act_psi_qiskit_string(circ, circuit.params, i_param, len(circuit.params))

        projector_circuits = []
        for j_param in range(len(hamiltonian.bonds)):
            circ_ham = circ.copy()
            circ_ham = circuit.act_permutation_qiskit(circ_ham, [hamiltonian.bonds[j_param]], ancilla=True, ancilla_qubit_idx = circuit.n_qubits)

            #print(hamiltonian.bonds[j_param])            
            for idx, cycl in enumerate(projector.cycl):
                if idx not in range(idx_from, idx_to):
                    continue

                current_circ = circ_ham.copy('circuit_projector_{:d}_{:d}'.format(j_param, idx))

                pair_permutations = projector.all_pair_permutations[idx]
                #print(pair_permutations)
                current_circ = circuit.act_permutation_qiskit(current_circ, pair_permutations, ancilla=True, ancilla_qubit_idx = circuit.n_qubits)
                current_circ.h(circuit.n_qubits)

                current_circ.measure(circuit.n_qubits, 0)
                projector_circuits.append(current_circ.copy())

        t = time()
        result = qiskit.execute(projector_circuits, backend, shots=N_samples).result()
        print('noisy simulation of der entry with quasm took', time() - t)

        survivals = []
        for i in range(len(projector_circuits)):
            survivals.append(-1. + 2. * result.data(i)['counts']['0x0'] / N_samples if '0x0' in result.data(i)['counts'].keys() else -1)

        survivals = np.array(survivals).reshape((len(hamiltonian.bonds), -1))
        survivals = np.einsum('ij,j->ij', survivals, np.array(projector.characters)).mean(axis=-1)
        der.append(survivals.dot(np.array(hamiltonian.js)))
        print(der)

    return -np.array(der)


def compute_energy_qiskit_hadamardtest(circuit, hamiltonian, projector, N_samples, noise_model):
    backend = qiskit.providers.aer.QasmSimulator(method='statevector', noise_model=noise_model, \
                                                 max_parallel_experiments=len(projector.cycl) * len(hamiltonian.bonds))

    circ = circuit.init_circuit_qiskit(ancilla=True)
    circ.h(circuit.n_qubits)
    circ = circuit.act_dimerization_qiskit(circ)
    circ = circuit.act_psi_qiskit_string(circ, circuit.params, 0, len(circuit.params))
    projector_circuits = []

    for j_param in range(len(hamiltonian.bonds)):
        circ_ham = circ.copy()
        circ_ham = circuit.act_permutation_qiskit(circ_ham, [hamiltonian.bonds[j_param]], ancilla=True, ancilla_qubit_idx = circuit.n_qubits)


        for idx, cycl in enumerate(projector.cycl):
            current_circ = circ_ham.copy('circuit_projector_{:d}_{:d}'.format(j_param, idx))

            pair_permutations = projector.all_pair_permutations[idx]
            current_circ = circuit.act_permutation_qiskit(current_circ, pair_permutations, ancilla=True, ancilla_qubit_idx = circuit.n_qubits)
            current_circ.h(circuit.n_qubits)

            current_circ.measure(circuit.n_qubits, 0)
            projector_circuits.append(current_circ.copy())

    t = time()
    result = qiskit.execute(projector_circuits, backend, shots=N_samples).result()
    print('noisy simulation of energy with quasm took', time() - t)

    survivals = []
    for i in range(len(projector_circuits)):
        survivals.append(-1. + 2. * result.data(i)['counts']['0x0'] / N_samples if '0x0' in result.data(i)['counts'].keys() else -1)

    survivals = np.array(survivals).reshape((len(hamiltonian.bonds), -1))
    survivals = np.einsum('ij,j->ij', survivals, np.array(projector.characters)).mean(axis=-1)

    return survivals.dot(np.array(hamiltonian.js))


def compute_gij_qiskit_hadamardtest(circuit, projector, N_samples, noise_model, real_imag, idx_from, idx_to):
    backend = qiskit.providers.aer.QasmSimulator(method='statevector', noise_model=noise_model, \
                                                 max_parallel_experiments=len(projector.cycl) * len(circuit.params) ** 2)

    gij = np.zeros((len(circuit.params), len(circuit.params)), dtype=np.complex128)
    t_total = time()
    projector_circuits = []
    for i_param in range(len(circuit.params)):
        for j_param in range(i_param, len(circuit.params)):
            circ = circuit.init_circuit_qiskit(ancilla=True)
            circ.h(circuit.n_qubits)
            circ = circuit.act_dimerization_qiskit(circ)
            circ = circuit.act_psi_qiskit_string(circ, circuit.params, 0, i_param)
            circ = circuit.act_derivative_qiskit(circ, i_param, ancilla_qubit_idx = circuit.n_qubits)
            circ.x(circuit.n_qubits)

            circ = circuit.act_psi_qiskit_string(circ, circuit.params, i_param, j_param)
            circ = circuit.act_derivative_qiskit(circ, j_param, ancilla_qubit_idx = circuit.n_qubits)
            circ = circuit.act_psi_qiskit_string(circ, circuit.params, j_param, len(circuit.params))
            circ.x(circuit.n_qubits)

            phase = 0 if real_imag == 'real' else +np.pi / 2.
            circ.u1(phase, circuit.n_qubits)


            for idx, cycl in enumerate(projector.cycl):
                if idx not in range(idx_from, idx_to):
                    continue
                current_circ = circ.copy('circuit_projector_{:d}'.format(idx))

                pair_permutations = projector.all_pair_permutations[idx]
                current_circ = circuit.act_permutation_qiskit(current_circ, pair_permutations, ancilla=True, ancilla_qubit_idx = circuit.n_qubits)
                current_circ.h(circuit.n_qubits)

                current_circ.measure(circuit.n_qubits, 0)
                projector_circuits.append(current_circ.copy())

    t = time()
    result = qiskit.execute(projector_circuits, backend, shots=N_samples).result()
    print('noisy simulation of all gij with quasm took', time() - t)

    survivals = []
    for i in range(len(projector_circuits)):
        survivals.append(-1. + 2. * result.data(i)['counts']['0x0'] / N_samples if '0x0' in result.data(i)['counts'].keys() else -1)

    ctr = 0
    for i_param in range(len(circuit.params)):
        for j_param in range(len(circuit.params)):
            for proj in range(idx_from, idx_to):
                if j_param < i_param:
                    continue
                gij[i_param, j_param] += survivals[ctr] * projector.characters[proj]
                ctr += 1
    gij /= len(range(idx_from, idx_to)) #[i_param, j_param] = np.mean(survivals)

    print('total time', time() - t_total)
    for i in range(len(circuit.params)):
        for j in range(i):
            if real_imag == 'real':
                gij[i, j] = gij[j, i]
            else:
                gij[i, j] = -gij[j, i]

    return gij

'''
def compute_energy_qiskit_hadamardtest(circuit, projector, hamiltonian, N_samples, n_repetitions, noise_model):
    circ = circuit.init_circuit_qiskit(ancilla=True)
    circ.h(circuit.n_qubits)
    circ = circuit.act_dimerization_qiskit(circ)
    circ = circuit.act_psi_qiskit(circ, circuit.params)


    survivals = np.zeros((len(projector.cycl), len(hamiltonian.bonds), n_repetitions), dtype=np.float64)
    for idx, cycl in enumerate(projector.cycl):
        projector_circuits = []
        current_circ = circ.copy('circuit_projector_{:d}'.format(idx))

        pair_permutations = projector.all_pair_permutations[idx]
        current_circ = circuit.act_permutation_qiskit(current_circ, pair_permutations, ancilla=True, ancilla_qubit_idx = circuit.n_qubits)

        for en_idx in range(len(hamiltonian.bonds)):
            name = 'circuit_projector_{:d}_energy_{:d}'.format(idx, en_idx)
            bond = hamiltonian.bonds[en_idx]
            j = hamiltonian.js[en_idx]

            current_circ_en = current_circ.copy(name)
            current_circ_en = circuit.act_permutation_qiskit(current_circ_en, [bond], ancilla=True, ancilla_qubit_idx = circuit.n_qubits)
            current_circ_en.h(circuit.n_qubits)
            current_circ_en.save_statevector(pershot=True)
            projector_circuits.append(current_circ_en.copy(name))

        t = time()

        backend = qiskit.providers.aer.QasmSimulator(method='statevector', noise_model=noise_model)
        result = qiskit.execute(projector_circuits, backend, shots=n_repetitions, memory=True, basis_gates = backend.configuration().basis_gates).result()
        print('noisy simulation with quasm took', time() - t)

        ctr = 0
        for j in range(len(hamiltonian.bonds)):
            survivals[idx, j, :] = np.sum(np.abs(np.array(result.data(ctr)['statevector'])[:, 2 ** (circuit.n_qubits):]) ** 2, axis=-1)
            ctr += 1

    
    survivals = np.clip(survivals, a_min=0., a_max=1.)
'''
def compute_gij_qiskit_survivalrate(circuit, new_parameters, N_samples, n_repetitions, noise_model):
    basis_gates = noise_model.basis_gates

    circ = circuit.init_circuit_qiskit()
    circ = circuit.act_dimerization_qiskit(circ)
    circ = circuit.act_psi_qiskit(circ, circuit.params)

    circ = circuit.act_psi_qiskit(circ, new_parameters, inverse=True)
    circ = circuit.act_dimerization_qiskit(circ, inverse=True)
    circ.save_statevector(pershot=True)

    t = time()

    backend = qiskit.providers.aer.QasmSimulator(method='statevector', noise_model=noise_model)
    result = qiskit.execute([circ], backend, shots=n_repetitions, basis_gates = backend.configuration().basis_gates).result()
    print('noisy simulation with quasm took', time() - t)

    survival = np.abs(np.array(result.data(0)['statevector'])[:, 0]) ** 2
    survival = np.clip(survival, a_min=0., a_max=1.)
    
    survival = np.random.binomial(N_samples, survival) / N_samples

    return np.mean(survival)




def compute_wavefunction(circuit, n_repetitions, noise_model):  # no projector allowed in this scenario
    circ = circuit.init_circuit_qiskit(ancilla = False)
    circ = circuit.act_dimerization_qiskit(circ)
    circ = circuit.act_psi_qiskit(circ, circuit.params)
    circ.save_statevector(pershot=True)


    backend = qiskit.providers.aer.QasmSimulator(method='statevector', noise_model=noise_model)
    result = qiskit.execute([circ], backend, shots=n_repetitions, memory=True, basis_gates = backend.configuration().basis_gates).result()

    wfs = np.array(result.data(0)['statevector']) #.mean(axis=0)
    return wfs


def compute_energy_qiskit_fromshift(circuit, hamiltonian, N_samples, n_repetitions, noise_model):  # no projector allowed in this scenario
    circ = circuit.init_circuit_qiskit(ancilla = False)
    circ = circuit.act_dimerization_qiskit(circ)
    circ = circuit.act_psi_qiskit(circ, circuit.params)
    circ.save_statevector(pershot=True)


    backend = qiskit.providers.aer.QasmSimulator(method='statevector', noise_model=noise_model)
    result = qiskit.execute([circ], backend, shots=n_repetitions, memory=True, basis_gates = backend.configuration().basis_gates).result()

    wfs = np.array(result.data(0)['statevector'])

    probas = np.zeros((n_repetitions, len(hamiltonian.js), 4), dtype=np.float64)

    uu = np.kron(np.diag([1, 0]), np.diag([1, 0]))
    ud = np.kron(np.diag([1, 0]), np.diag([0, 1]))
    du = np.kron(np.diag([0, 1]), np.diag([1, 0]))
    dd = np.kron(np.diag([0, 1]), np.diag([0, 1]))

    for idx_op, op in enumerate([uu, ud, du, dd]):
        for idx_term, term in enumerate(hamiltonian.bonds):
            OP = ls.Operator(circuit.basis, [ls.Interaction(op, [term])])
            probas[:, idx_term, idx_op] = np.einsum('ji,ji->i', wfs.T.conj(), OP(wfs.T)).real

    probas_sampled = np.random.binomial(N_samples, probas) / N_samples
    probas_sampled = probas_sampled[..., 0] + probas_sampled[..., 3] - probas_sampled[..., 1] - probas_sampled[..., 2]
    probas_sampled = probas_sampled.mean(axis=0) * 3  # SU(2) invariance
    return probas_sampled.dot(np.array(hamiltonian.js))



