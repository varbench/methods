import numpy as np
from scipy.optimize import minimize
import qutip as qu
from qutip.qobj import Qobj
from qutip.operators import identity
from qutip.qip.operations import gate_expand_2toN, rx, ry,cnot
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

print("Heuristic RY-CNOT circuit")

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
    psi_initial = qu.basis(2,1)
    for i in range(n-1):
        psi_initial = qu.tensor(psi_initial, qu.basis(2,1))

    return psi_initial/psi_initial.norm()

def blockCNOT(n_qubits,pairs):
    '''a ladder of CNOTS gates'''
    entangler = 1

    operator_full = qu.Qobj()

    for term in range(len(pairs)):

         entangler= cnot(n_qubits,control=pairs[term][0],target=pairs[term][1])*entangler

    return entangler


def rotation_blockY(n_qubits, thetas,list):
    '''a block of parametrized RY gates'''

    thetas_counter = 0


    if 0 in list :
        z_block = ry(thetas[thetas_counter])
        thetas_counter += 1
    else:
        z_block = qu.qeye(2)

    for i in range(1, n_qubits):
        if i in list :
            z_block = tensor(z_block,ry(thetas[thetas_counter]))
            thetas_counter += 1
        else:
            z_block = tensor(z_block,qu.qeye(2))


    return z_block

def ansatzRYCNOT(thetas, n_qubits, hamiltonian, psi_in):
    '''the circuit ansatz'''

    psi_circuit = qu.identity(2)
    for i in range(1, n_qubits):
        psi_circuit = tensor(psi_circuit, identity(2))

    thetas_counter = 0

    #defining the qubits on which the RY gates and the CNOTs act
    pairs = []
    list=[]
    for i in range(n_qubits-1):
        pairs.append([i,i+1])
    for i in range(n_qubits):
        list.append(n_qubits-i-1)


    # create the ansatz with the custom layers defined above
    for block in range(n_blocks):

        if( block > 0):
            psi_in=psi_circuit


        psi_circuit = rotation_blockY(n_qubits=n_qubits,thetas=thetas[thetas_counter:thetas_counter + len(list)],list=list)*psi_in
        thetas_counter += len(list)

        psi_circuit = blockCNOT(n_qubits=n_qubits,pairs=pairs)*psi_circuit

    psi_circuit = rotation_blockY(n_qubits=n_qubits,thetas=thetas[thetas_counter:thetas_counter +len(list)],list=list)*psi_circuit

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
n_blocks = 10
print("Heuristic RY-CNOT ansatz with", n_blocks, "layers")


psi_in = initstate(n_qubits)


print("Loading the WF params\n")
# hardcoding the best parameters found for this ansatz (with 10 layers)
thetas=[ 4.15245693e-01, -6.92705758e-01, -8.01379784e-01, -1.04233504e+00,
       -6.66588832e-01,  1.01800203e+00,  1.56600807e+00,  1.55997659e+00,
        1.57516196e+00, -1.56922831e+00,  2.53732852e+00, -1.57151849e+00,
       -6.81159310e-01,  2.17771831e+00, -9.28620757e-01, -1.59808451e+00,
        3.12743849e+00, -4.46739672e-02,  8.17467775e-04,  1.32057389e-03,
        2.68895995e+00, -2.70637649e+00,  2.05823506e-01, -4.27482382e+00,
        1.63964612e+00, -1.62843431e+00,  8.58274213e-02,  3.16124363e+00,
       -3.00142053e-03, -1.17813578e-03,  3.14188290e+00,  3.25677350e-02,
       -4.08508660e+00, -1.12990374e+00, -4.66085521e-01,  2.47257964e-02,
        3.03924565e+00,  3.16605936e+00,  3.13864120e+00, -1.98953910e-03,
        3.91490701e-06, -2.24074124e-01, -1.77681918e+00, -1.55350387e+00,
       -1.64938759e+00,  1.09150494e+00, -2.48374492e+00, -1.12611605e-02,
       -1.25347252e-02,  1.84878526e-04, -3.14172575e+00,  2.68050823e+00,
       -2.23803092e+00, -1.19588468e+00,  2.49867559e-01, -3.30028453e+00,
       -1.64435928e+00,  8.03282582e-02, -9.02800927e-02,  7.71669075e-05,
       -5.75740621e-05, -2.77150272e+00, -2.33207844e+00, -2.39547544e+00,
       -1.54659784e+00, -1.55526274e+00, -1.54152864e+00,  2.25015253e-01,
       -1.53882112e-02, -3.79646606e-02,  3.14153827e+00, -1.64344109e-02,
       -5.33756473e-02,  4.35241913e-01,  9.13860730e-02,  2.31558754e+00,
        1.56576474e+00, -6.09712908e-02,  3.11565663e-01, -1.25115642e-02,
        2.59911626e-05,  1.12092294e-02,  1.04005108e-01,  4.37874047e-02,
       -3.29435192e-01, -2.57317596e-01, -6.08383485e-02, -1.02529075e+00,
        1.77635511e-01, -1.25724666e-01,  1.57077676e+00,  1.05154537e+00,
        1.15057226e+00, -1.84657436e+00, -1.03619320e+00, -8.34276657e-01,
        7.56500968e-01,  2.94362069e-01, -9.17551484e-01, -5.56706030e-01,
        2.78196401e-05,  6.03530889e-05,  5.59487193e-05,  4.08932190e-05,
        3.35448107e-05,  2.04164586e-05,  1.81066126e-05, -2.80071421e-05,
       -1.55397564e-05, -1.37363070e-05]



# evaluate the expectation values of H and H^2 on the loaded ansatz
val, var = ansatzRYCNOT(thetas, n_qubits, hamiltonian, psi_in)
print("Energy", val, ", Variance", var)

vscore = n_qubits*var / (val)**2
rel_err = np.abs((val - exacten)/exacten)
print("Vscore", vscore, ", Relative energy error", rel_err)
