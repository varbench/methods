import numpy as np
from scipy.optimize import minimize
import qutip as qu
from qutip.qobj import Qobj
from qutip.operators import identity
from qutip.qip.operations import gate_expand_2toN, rx
from qutip.tensor import tensor
import sys

__author__ = 'Guglielmo Mazzola'


'''
Works with qutip 4.7.0
'''
print("Your version of qutip is", qu.__version__)

##################################

'''***System parameters***'''
n=10 #number of sites
gamma=-1.0
J=-1.0
pbc=False #open boundary conditions

print("Hamiltonian variational circuit")

###############################
'''***Functions***'''

def hamiltonian_TF(n,gamma,J,pbc):
    '''create a TF Ising hamiltonian, used to compute expectation values and exact gs'''

    dims = [[2] * n, [2] * n]
    pairs = []

    for i in range(n-1):
            pairs.append([i,i+1])

    if(pbc): pairs.append([0,i+1])

    operator_full = Qobj()

    for term in range(len(pairs)):

        if 0 in pairs[term]:
            operator_ij = qu.sigmaz()
        else:
            operator_ij = qu.qeye(2)

        for nn in range(1, n):

            if nn in pairs[term]:

                operator_ij = qu.tensor(operator_ij,  qu.sigmaz())

            else:
                operator_ij = qu.tensor(operator_ij, qu.qeye(2))

        operator_full +=   J*operator_ij

    for iop in range(n):

        if iop==0:
            operator_ij = qu.sigmax()
        else:
            operator_ij = qu.qeye(2)

        for nn in range(1, n):
            if iop==nn:
                operator_ij = qu.tensor(operator_ij,  qu.sigmax())
            else:
                operator_ij = qu.tensor(operator_ij, qu.qeye(2))

        operator_full +=   gamma*operator_ij

    return Qobj(operator_full, dims=dims)




def initstate(L):
    '''initialize in superposition of single qubit |+> states'''
    psi_initial = qu.basis(2,0)+qu.basis(2,1)
    for i in range(L-1):
        psi_initial = tensor(psi_initial, qu.basis(2,0)+qu.basis(2,1))
    return psi_initial/psi_initial.norm()



def myZZ(arg_value, N=None, control=0, target=1):
    '''custom exp(i theta ZZ) two qubit gate'''
    if (control == 1 and target == 0) and N is None:
        N = 2

    if N is not None:
        return gate_expand_2toN(myZZ(arg_value), N, control, target)

    return Qobj([[ np.exp(arg_value*1j), 0., 0,  0.],
                        [0, np.exp(-arg_value*1j),  0., 0],
                        [0.,  0., np.exp(-arg_value*1j), 0.],
                        [ 0., 0, 0,  np.exp(arg_value*1j)]],
                        dims=[[2, 2], [2, 2]])


def entangler_blockZZ(n_qubits, theta):
    '''one block of ZZ gates (with same theta) over the qubit register'''
    entangler = 1
    for term in range(n_qubits-1):
         entangler= myZZ(theta,n_qubits,control=term,target=term+1)*entangler

    return entangler


def rotation_blockX(n_qubits, theta):
    '''one block of RX gates (with same theta) over the qubit register'''
    x_block = rx(theta)
    for i in range(1, n_qubits):
        x_block = tensor(x_block,rx(theta))

    return x_block


def ansatzHV(thetas, n_qubits, hamiltonian, psi_in):
    '''the circuit ansatz'''

    psi_circuit = qu.identity(2)
    for i in range(1, n_qubits):
        psi_circuit = tensor(psi_circuit, identity(2))

    thetas_counter = 0
    # create the ansatz with the custom layers defined above
    for block in range(n_blocks):
        if( block == 0):
            psi_in=psi_in
        if( block > 0):
            psi_in=psi_circuit

        psi_circuit = entangler_blockZZ(n_qubits,thetas[thetas_counter])*psi_in
        thetas_counter += 1
        psi_circuit = rotation_blockX(n_qubits,thetas[thetas_counter])*psi_circuit
        thetas_counter += 1

    # measure the expectation value
    val = qu.expect(hamiltonian, psi_circuit)

    Hsquared = hamiltonian*hamiltonian
    H2expect = qu.expect(Hsquared,psi_circuit)


    var=H2expect- val*val


    return val, var

############################################

'''Exact diagonalization'''

hamiltonian = hamiltonian_TF(n,gamma,J,pbc)
n_qubits = len(hamiltonian.dims[0])
exacten = hamiltonian.eigenstates()[0][0]
print("Exact energy", exacten, ", L =",n_qubits)



'''***Circuit parameters***'''
n_blocks = 24
print("Hamiltonian variational ansatz with", n_blocks, "layers")


psi_in = initstate(n_qubits)

print("Loading the WF params\n")
# hardcoding the best parameters found for this ansatz (with 24 layers)
thetas=np.zeros(48)
thetas[:48]=[ 1.38944938, -2.68598443, -0.89050956,  0.74384366, -0.58385548,
        0.21147651,  0.79150307, -0.58287565, -0.54710389,  0.83314987,
        1.22132212, -2.61867202,  0.71995241,  0.3716459 , -1.66154474,
        2.1720757 ,  0.58037918, -1.83128399,  0.90047899, -1.97616278,
        0.29164909, -1.58402846,  0.36512155, -2.94838537, -0.23599522,
        2.12423786, -0.26066222,  1.46717049, -0.87275616,  1.94350639,
       -0.71357013,  2.6799029 , -0.53164801,  0.27716476,  0.47967709,
        0.35666815,  0.43110715, -0.14770922, -0.26459335, -2.00376091,
       -0.14410481,  1.2425756 ,  0.32887615, -0.37830254,  0.06504382,
       -0.54147089,  0.01415658,  0.00697829]



# evaluate the expectation values of H and H^2 on the loaded ansatz
val, var = ansatzHV(thetas, n_qubits, hamiltonian, psi_in)
print("Energy", val, ", Variance", var)

vscore = n_qubits*var / (val)**2
rel_err = np.abs((val - exacten)/exacten)
print("Vscore", vscore, ", Relative energy error", rel_err)
