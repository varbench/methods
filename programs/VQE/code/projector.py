import numpy as np 
import lattice_symmetries as ls
import utils
from copy import deepcopy

class Projector(object):
    def __call__(self, state, n_term = None, inv=False):
        if n_term is not None:
            if inv:
                return state[..., self.permutations_inv[n_term]] * 1. / self.characters[n_term] if self.characters[n_term] != 1.0 else state[..., self.permutations_inv[n_term]]
            return state.T[self.permutations[n_term], ...].T * self.characters[n_term] if self.characters[n_term] != 1.0 else np.ascontiguousarray(state.T)[self.permutations[n_term], ...].T

        state_projected = state * 0.0
        if not inv:
            for permutation, character in zip(self.permutations, self.characters):
                state_projected += state[..., permutation] * character
        else:
            for permutation, character in zip(self.permutations_inv, self.characters):
                state_projected += state[..., permutation] * 1. / character
        return state_projected / len(self.permutations)
        

class ProjectorFull(Projector):
    def __init__(self, Hbonds, n_qubits, su2, basis, generators, eigenvalues, degrees):
        self.basis = basis
        self.n_qubits = n_qubits
        self.basis_size = basis.number_states
        self.Hbonds = Hbonds

        self.maps, self.permutations, self.characters, self.cycl, self.all_pair_permutations, self.all_double_permutations = self._init_projector(generators, eigenvalues, degrees)
        self.permutations_inv = [np.argsort(perm) for perm in self.permutations]
        self.nterms = len(self.permutations)

        return


    def _combine_generators(self, generators, eigenvalues, degrees):
        permutations = [np.arange(self.basis_size)]
        maps = [np.arange(self.n_qubits)]
        characters = [1. + 0.0j]

        if len(generators) == 0:
            return permutations, maps, characters

        for generator, eigenvalue, degree in zip(generators, eigenvalues, degrees):
            m, g = generator
            m = m.copy()
            g = g.copy()

            lamb = eigenvalue

            new_permutations = []
            new_characters = []
            new_maps = []

            for d in range(degree - 1):
                permutations_d = []
                characters_d = []
                maps_d = []
                for mm, symm, ch in zip(maps, permutations, characters):
                    permutations_d.append(symm[g])
                    characters_d.append(lamb * ch)
                    maps_d.append(mm[m])

                new_permutations.append(deepcopy(permutations_d))
                new_characters.append(deepcopy(characters_d))
                new_maps.append(deepcopy(maps_d))

                lamb *= eigenvalue
                g = g[generator[1]]
                m = m[generator[0]]

            for perm in new_permutations:
                permutations += perm
            for ch in new_characters:
                characters += ch
            for x in new_maps:
                maps += x
        return permutations, maps, characters

    def _init_projector(self, generators, eigenvalues, degrees):
        permutations, maps, characters = self._combine_generators(generators, eigenvalues, degrees)

        total = 1
        for d in degrees:
            total *= d


        assert len(permutations) == total

        print('terms in the projector:', total)


        self.list_of_pairs = []
        added = np.zeros(len(maps))
        friendless = 0
        for idx, p in enumerate(maps):
            friend = False
            if np.allclose(np.arange(len(maps[0])), p[p]):
                friendless += 1
            for k in range(idx + 1, len(maps)):
                if np.allclose(p, maps[k]):
                    print('coincide', idx, k, p)
                    total -= 1
                if np.allclose(p, np.argsort(maps[k])):
                    print(idx, k, 'are friends')
                    self.list_of_pairs.append((idx, k))
                    added[idx] = 1
                    added[k] = 1
                    friend = True
            if not friend:
                print(idx, 'is friendless')
                if not added[idx]:
                    self.list_of_pairs.append((idx,))

        print(self.list_of_pairs)
        print('unique', total)
        assert sum([len(x) for x in self.list_of_pairs]) == len(maps)
        print('friendless', friendless)


        ### speed-up for s3it ###
        if len(generators) == 4:
            self.lpermutations, self.lmaps, self.lcharacters = self._combine_generators(generators[:2], eigenvalues[:2], degrees[:2]) 
            self.rpermutations, self.rmaps, self.rcharacters = self._combine_generators(generators[2:], eigenvalues[2:], degrees[2:])
        elif len(generators) == 3:
            self.lpermutations, self.lmaps, self.lcharacters = self._combine_generators(generators[:1], eigenvalues[:1], degrees[:1])
            self.rpermutations, self.rmaps, self.rcharacters = self._combine_generators(generators[1:], eigenvalues[1:], degrees[1:])
        else:
            self.lpermutations, self.lmaps, self.lcharacters = self._combine_generators(generators, eigenvalues, degrees)
            self.rpermutations, self.rmaps, self.rcharacters = self._combine_generators(generators, eigenvalues, degrees)


        self.left_right_decompositions = []	
        for idx, p in enumerate(maps):
            found = False
            for idxl, l in enumerate(self.lmaps):
                for idxr, r in enumerate(self.rmaps):
                    if np.allclose(p, np.argsort(l)[r]) and not found:
                        self.left_right_decompositions.append((idxl, idxr))
                        found = True
        print(self.left_right_decompositions)
        assert len(self.left_right_decompositions) == len(maps)


        def cycles(perm):
            remain = set(perm)
            result = []
            while len(remain) > 0:
                n = remain.pop()
                cycle = [n]
                while True:
                    n = perm[n]
                    if n not in remain:
                        break
                    remain.remove(n)
                    cycle.append(n)
                result.append(cycle)
            return result


        cycl = []
        for idx, p in enumerate(maps):
            map_swaps = []
            cycl.append(cycles(p))


        all_pair_permutations = []
        for c in cycl:
            pair_permutations = []

            for loop in c:
                if len(loop) == 1:
                    continue

                for i in reversed(range(1, len(loop))):
                    pair_permutations.append((loop[0], loop[i]))
            all_pair_permutations.append(deepcopy(pair_permutations))



        all_double_permutations = []
        '''
        for pair_permutations, permutation, m in zip(all_pair_permutations, permutations, maps):
            double_permutations = []
            for separator in range(len(pair_permutations)):
                perm_before = np.arange(self.n_qubits)
                perm_after = np.arange(self.n_qubits)

                for pair in pair_permutations[:separator]:
                    perm_before[pair[0]], perm_before[pair[1]] = perm_before[pair[1]], perm_before[pair[0]]
                for pair in pair_permutations[separator:]:
                    perm_after[pair[0]], perm_after[pair[1]] = perm_after[pair[1]], perm_after[pair[0]]
                #perm_before = np.argsort(perm_before)
                #perm_after = np.argsort(perm_after)

                permutation_before = np.argsort(utils.spin_to_index(utils.index_to_spin(self.basis.states, number_spins = self.n_qubits)[:, perm_before], number_spins = self.n_qubits)) 
                permutation_after = np.argsort(utils.spin_to_index(utils.index_to_spin(self.basis.states, number_spins = self.n_qubits)[:, perm_after], number_spins = self.n_qubits))
                double_permutations.append((permutation_before.copy(), permutation_after.copy()))

                assert np.allclose(perm_before[perm_after], m)
                assert np.allclose(permutation, permutation_before[permutation_after])
            all_double_permutations.append(deepcopy(double_permutations))
        '''



        return maps, permutations, characters, cycl, all_pair_permutations, all_double_permutations



