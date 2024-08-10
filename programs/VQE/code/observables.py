import numpy as np
import lattice_symmetries as ls
import os



sz = np.array([[1, 0], \
               [0, -1]])

sx = np.array(
               [[0, 1], \
               [1, 0]]
              )

sy = np.array([[0, 1.0j], [-1.0j, 0]])

SS = np.kron(sx, sx) + np.kron(sy, sy) + np.kron(sz, sz)

def neel_order(Lx, Ly, basis, su2=False):
    n_qubits = Lx * Ly

    site_plus = []
    site_minus = []
    for x in range(Lx):
        for y in range(Ly):
            if (x + y) % 2 == 0:
                site_plus.append(x + y * Lx)
            else:
                site_minus.append(x + y * Lx)

    return ls.Operator(basis, [ls.Interaction(sz / n_qubits, site_plus), ls.Interaction(-sz / n_qubits, site_minus)]), 'Neel'

def neel_order_hexagon(basis, su2=False):
    n_qubits = 6

    site_plus = []
    site_minus = []
    for i in range(n_qubits):
        if i % 2 == 0:
            site_plus.append(i)
        else:
            site_minus.append(i)

    return ls.Operator(basis, [ls.Interaction(sz / n_qubits, site_plus), ls.Interaction(-sz / n_qubits, site_minus)]), 'Neel_hex'


def neel_order_honeycomb(Lx, Ly, basis, su2=False):
    n_qubits = 2 * Lx * Ly

    site_plus = []
    site_minus = []
    for i in range(n_qubits):
        if i % 2 == 0:
            site_plus.append(i)
        else:
            site_minus.append(i)

    return ls.Operator(basis, [ls.Interaction(sz / n_qubits, site_plus), ls.Interaction(-sz / n_qubits, site_minus)]), 'Neel_honey'



def stripe_order(Lx, Ly, basis, su2=False):
    n_qubits = Lx * Ly

    site_plus = []
    site_minus = []
    for x in range(Lx):
        for y in range(Ly):
            if x % 2 == 0:
                site_plus.append(x + y * Lx)
            else:
                site_minus.append(x + y * Lx)

    return ls.Operator(basis, [ls.Interaction(sz / n_qubits, site_plus), ls.Interaction(-sz / n_qubits, site_minus)]), 'Stripe'


def dimer_order(Lx, Ly, basis, su2=False, BC='PBC'):
    n_qubits = Lx * Ly
    

    bond_plus = []
    bond_minus = []
    for x in range(Lx):
        for y in range(Ly):
            if x % 2 == 0:
                if x < Lx - 1 or BC == 'PBC':
                    bond_plus.append(((x + y * Lx), (((x + 1) % Lx) + y * Lx)))
            else:
                if x < Lx - 1 or BC == 'PBC':
                    bond_minus.append(((x + y * Lx), (((x + 1) % Lx) + y * Lx)))

    return ls.Operator(basis, [ls.Interaction(SS / n_qubits, bond_plus), ls.Interaction(-SS / n_qubits, bond_minus)]), 'Dimer'



def dimer_order_honeycomb(Lx, Ly, basis, su2=False, BC='PBC'):
    n_qubits = 2 * Lx * Ly

    bonds = []
    for i in range(n_qubits // 2):
        bonds.append((2 * i, 2 * i + 1))
    return ls.Operator(basis, [ls.Interaction(SS / n_qubits, bonds)]), 'dimer_honeycomb'



def plaquette_order_hexagon(basis, su2=False, BC='PBC'):
    n_qubits = 6

    bond_plus = [(0, 1), (2, 3), (4, 5)]
    bond_minus = [(1, 2), (3, 4), (0, 5)]

    return ls.Operator(basis, [ls.Interaction(SS / n_qubits, bond_plus), ls.Interaction(-SS / n_qubits, bond_minus)]), 'Plaquette_hex'


def plaquette_order_honeycomb_3x3(basis, su2=False, BC='PBC'):
    n_qubits = 18

    bond_plus = [(0, 1), (3, 6), (2, 15), (8, 9), (11, 14), (5, 10), (16, 17), (4, 13), (7, 12)]
    bond_minus = [(1, 6), (2, 3), (0, 15), (9, 14), (10, 11), (5, 8), (17, 4), (12, 13), (7, 16)]

    return ls.Operator(basis, [ls.Interaction(SS / n_qubits, bond_plus), ls.Interaction(-SS / n_qubits, bond_minus)]), 'Plaquette_hex_3x3'


class Observables(object):
    def __init__(self, config, hamiltonian, circuit, projector):
        self.path_to_logs = config.path_to_logs
        openmode = 'a' if config.mode == 'continue' else 'w'

        self.main_log = open(os.path.join(self.path_to_logs, 'main_log.dat'), openmode)
        self.force_log = open(os.path.join(self.path_to_logs, 'force_log.dat'), openmode)
        self.exact_force_log = open(os.path.join(self.path_to_logs, 'exact_force_log.dat'), openmode)
        self.force_SR_log = open(os.path.join(self.path_to_logs, 'force_SR_log.dat'), openmode)
        self.exact_force_SR_log = open(os.path.join(self.path_to_logs, 'exact_force_SR_log.dat'), openmode)

        self.parameters_log = open(os.path.join(self.path_to_logs, 'parameters_log.dat'), openmode)
        self.lambda_log = open(os.path.join(self.path_to_logs, 'lambda_log.dat'), openmode)

        self.observables = config.observables
        self.hamiltonian = hamiltonian
        self.projector = projector
        self.circuit = circuit
        self.config = config


        if config.write_logs:
            ### prepare main log ###
            string = 'gsenergy energy fidelity norm '
            for _, name in self.observables:
                string += name + ' '
            self.main_log.write(string + '\n')

        return

    def write_logs(self):
        #print('smth entered', flush=True)
        if self.config.test or self.config.N_samples is None:
            force_exact = self.circuit.forces_exact
            for f in force_exact:
                self.exact_force_log.write('{:.4f} '.format(f))
            self.exact_force_log.write('\n')
            self.exact_force_log.flush()

            force_SR_exact = self.circuit.forces_SR_exact
            for f in force_SR_exact:
                self.exact_force_SR_log.write('{:.4f} '.format(f))
            self.exact_force_SR_log.write('\n')
            self.exact_force_SR_log.flush()

        force = self.circuit.forces
        if force is not None:
            for f in force:
                self.force_log.write('{:.4f} '.format(f))
            self.force_log.write('\n')
            self.force_log.flush()


        force_SR = self.circuit.forces_SR
        if force_SR is not None:
            for f in force_SR:
                self.force_SR_log.write('{:.4f} '.format(f))
            self.force_SR_log.write('\n')
            self.force_SR_log.flush()


        parameters = self.circuit.params
        for p in parameters:
            self.parameters_log.write('{:.4f}, '.format(p))
        self.parameters_log.write('\n')
        self.parameters_log.flush()

        lamb = self.circuit.lamb
        self.lambda_log.write('{:.4f}\n'.format(lamb))
        self.lambda_log.flush()


        #### compute energy ###
        state = self.circuit(noisy=False)
        state_proj = self.projector(state)
        norm = np.vdot(state, state_proj)
        energy = (np.vdot(state, self.hamiltonian(state_proj)) / norm).real - self.hamiltonian.energy_renorm
        #print(np.vdot(state, self.hamiltonian(state_proj)) / norm)
        #print(norm, 'norm full')
        #print(np.vdot(self.projector(self.hamiltonian(state)), self.hamiltonian(self.projector(state))) / np.vdot(self.projector(self.hamiltonian(state)), self.projector(self.hamiltonian(state))))

        ### compute fidelity ###
        state_proj = state_proj / np.sqrt(np.vdot(state_proj, state_proj))
        #print(norm, np.vdot(state, state), np.vdot(self.hamiltonian.ground_state[0], self.hamiltonian.ground_state[0]))
        # exit(-1)
        assert np.isclose(np.vdot(state_proj, state_proj), 1.0)
        fidelity = np.abs(np.vdot(self.hamiltonian.ground_state[0], state_proj)) ** 2
        #for i, j in zip(self.hamiltonian.ground_state[0], state_proj):
        #    print(i, j)
        #print(fidelity)


        obs_vals = []
        obs_vals_ed = []
        for operator, _ in self.observables:
            val = np.dot(state_proj.conj(), operator(operator(state_proj)))
            val_ed = np.dot(self.hamiltonian.ground_state[0].conj(), operator(operator(self.hamiltonian.ground_state[0])))
            assert np.isclose(val.imag, 0.0)
            assert np.isclose(val_ed.imag, 0.0)
            obs_vals.append(val.real)
            obs_vals_ed.append(val_ed.real)
        
        #### dimer amendment ###
        dimer_avg = np.dot(state_proj.conj(), self.observables[-1][0](state_proj))
        assert np.isclose(dimer_avg.imag, 0.0)
        obs_vals[-1] -= (dimer_avg ** 2).real

        dimer_avg = np.dot(self.hamiltonian.ground_state[0].conj(), self.observables[-1][0](self.hamiltonian.ground_state[0]))
        assert np.isclose(dimer_avg.imag, 0.0)
        obs_vals_ed[-1] -= (dimer_avg ** 2).real

        data = [j for i in zip(obs_vals, obs_vals_ed) for j in i]

        self.main_log.write(('{:.7f} {:.7f} {:.14f} {:.7f} ' + '{:.7f}/{:.7f} ' * len(obs_vals) + '\n').format(self.hamiltonian.gse - self.hamiltonian.energy_renorm, energy, fidelity, norm.real, *data))
        print(('{:.7f} {:.7f} {:.14f} {:.7f} ' + '{:.7f}/{:.7f} ' * len(obs_vals) + '\n').format(self.hamiltonian.gse - self.hamiltonian.energy_renorm, energy, fidelity, norm.real, *data), flush=True)
        self.main_log.flush()
        #exit(-1)
