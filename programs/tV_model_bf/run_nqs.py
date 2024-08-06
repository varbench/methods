import numpy as np
import netket as nk
from netket.utils import mpi, HashableArray
import jax
import flax
import jax.numpy as jnp
from utils import get_max_conn
from other_obs import structure_factor, dens_dens_corr, renormalized_corr, corr_fourier

from utils import distances_rij, k_vector, k_orbitals

import sys
#sys.dont_write_bytecode = True

import os

# detect MPI rank
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
MPI4JAX_NO_WARN_JAX_VERSION=1

# set only one visible device
os.environ["CUDA_VISIBLE_DEVICES"] = f"{rank}"
# force to use gpu
os.environ["JAX_PLATFORM_NAME"] = "gpu"

print('rank, jax.devices')
print(f"{rank} -> {jax.devices()}")

from jax.config import config
config.update('jax_disable_jit', False)

from hamiltonian import t_v_model
#from models import MeanBackflowSlater
from symm_models import SymmRBM, SymmMeanBackflowSlater
from symmetries import _character


nk.config.netket_experimental_fft_autocorrelation = True

import argparse

# Create the argument parser
parser = argparse.ArgumentParser()

# Add positional arguments
parser.add_argument("--L", type=int, help="Lattice size", required=True)
parser.add_argument("--Nf", type=int, help="Number of fermions", required=True)
parser.add_argument("--V", type=float, help="Coulomb interaction", required=True)
parser.add_argument("--j", type=int, help="Jastrow", required=True)
#parser.add_argument("--symm", type=int, help="Symmetries", default=True)


#if parser.parse_known_args()[0].J:
    #parser.add_argument("--J_rbm", type=int, help="Jastrow_RBM", default=True)

# Add optional arguments for symm=True case
#if parser.parse_known_args()[0].symm:
    #parser.add_argument("--charac", type=int, help="Character", required=True)

parser.add_argument("--bf", type=int, help="NN Backflow")

if parser.parse_known_args()[0].bf:
    parser.add_argument("--depth", type=int, help="Depth", required=True)
    parser.add_argument("--feat", type=int, help="Features", required=True)
else:
    parser.add_argument("--depth", type=int, help="Depth", default=1)
    parser.add_argument("--feat", type=int, help="Features",default=1 )

# Parse the arguments
args = parser.parse_args()

# Access the values of the arguments
print("L:", args.L)
print("Nf:", args.Nf)
print("V:", args.V)
#print("symmetrize:", bool(args.symm))

print("Jastrow:", bool(args.j))



#if args.J:
    #print("Jastrow RBM:", bool(args.J_rbm))
#if args.symm:
    #print("character:", args.charac)
if args.bf: 
    print("backflow:", bool(args.bf))
    print("depth:", args.depth)
    print("features:", args.feat)

print()
print()


print("number of tasks = ", mpi.n_nodes)
##################################################################################################################

###### system model parameters ######

D = 2  # dimension D
L = args.L  # L^D lattice
t = 1. # hopping
V = args.V # Coulomb


Ns = args.L**D # number of sites
Nf = args.Nf # number of fermions
print('Number of fermions = ',Nf, '\n')



###### Hamiltonian construction ######
# construct t-V model with corresponding hilbert space, graph and hamiltonian
ferm_hub = t_v_model(L,D,t,V,Nf)
hi = ferm_hub[0]
g = ferm_hub[1]
ham = ferm_hub[2]
print('ham max conn size', ham.max_conn_size)

ham._max_conn_size=get_max_conn(ham, change=True)#ham.max_conn_size
#ham._setup()
print('new ham max conn size', ham._max_conn_size)
print('double check ham max conn size', ham.max_conn_size)

#print('get max conn function', get_max_conn(ham, change=True))

n_chains = 288 #(n_tasks)
n_samples=1024
chunk_size=None
n_sweeps=hi.size
learning_rate=0.025
diag_shift=0.01
n_samples_v = 4196*2
iterations = 1500


print('n_chains', n_chains)
print('n_samples', n_samples)
print('n_sweeps', n_sweeps)
print('learning_rate', learning_rate)
print('diag_shift', diag_shift)
print('n_samples_v', n_samples_v)



##################################################################################################################

def run_ansatz(i):

    # choose sampler
    sa = nk.sampler.MetropolisExchange(hi, graph=g, n_chains=n_chains,n_sweeps=n_sweeps)
    

    # define the model
    #if args.symm == 0:
        #print('no sym')
        #---- without lattice symmetries ----
        #ma = MeanBackflowSlater(L=args.L, D=D, Nf=args.Nf, Ns=Ns, mf_orbitals=True, backflow=args.bf, jastrow=args.J)
    
    #elif args.symm:
    #---- with lattice symmetries ----
    ma = SymmMeanBackflowSlater(L=args.L,D=D, Nf=args.Nf, Ns=Ns, symmetries = symm, character=character, graph=g,
                                backflow=args.bf, mf_orbitals=True, jastrow=args.j, depth=args.depth, features = args.feat)
    
    # variational state
    vs = nk.vqs.MCState(sa, ma, n_samples=n_samples,chunk_size=chunk_size, n_discard_per_chain=2)#chunk_size=8)#n_discard_per_chain=100
    vs.init_parameters()
    print('acceptance', vs.sampler_state.n_accepted_proc, vs.sampler_state.acceptance)
    

    print('number of params in the model', vs.n_parameters,)
    #np.savetxt(f"output/n_params.txt",[vs.n_parameters])

    # print model architectzre and the number of samples per rank
    print(jax.tree_util.tree_map(lambda x: f"{x.shape} -> {x.dtype}", vs.parameters))
    print('samples per rank', vs.n_samples_per_rank)


    # use sgd with Stochastic Reconfiguration
    opt = nk.optimizer.Sgd(learning_rate=learning_rate)
    qgt = nk.optimizer.qgt.QGTJacobianDense(holomorphic=False)#QGTOnTheFly()
    sr = nk.optimizer.SR(qgt, diag_shift=diag_shift)

    # VMC driver
    gs = nk.driver.VMC(ham, opt, variational_state=vs, preconditioner=sr)
    # run the optimization
    out_name = f"output/out_{i}"
    gs.run(n_iter=iterations, out=out_name, callback=acceptance_callback, )#obs={'Structure Factor': s})

    

    # save variational state prams --> flax serialization
    with open(f"output/out_{i}_params.mpack", 'wb') as file:
        file.write(flax.serialization.to_bytes(vs.parameters))
    

    return gs.state


def mean_en_var(ham, vs):

    print('check again ham', ham.max_conn_size)
    
    n_samples = n_samples_v # n_samples needs to be divisble by chunk_size (and number of samples per rank)

    vs.n_samples = n_samples
    vs.chunk_size=4196
    energy = vs.expect(ham)
    energy_mean = energy.mean.real

    variance = energy.variance.real
    #var_tot=s2_tot-en_tot**2
    
    #s = sqrt(var_tot / float(n_samples))
    standard = energy.error_of_mean
    #s2_tot = energy.variance.real + energy.mean.real ** 2.
    
    print('Energy, sigma, var')
    print(energy_mean, standard, variance)


    
    


    ####### k values (1D) ######
    k =  2 * np.pi / (L + 1) * np.arange(0, L + 1) # same as numpy fft
    #k = np.array([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    #print('k', k)
    ####### K values (2D) ######
    K = k_vector(L)
    print('Ks', np.array(K))

    
    ##############################################   
    ###### density-density correlation function ######
    print()
    C_r = []
    for r in range(0,L+1):
        corr = dens_dens_corr(g, Nf, r,)
        corr =  vs.expect(corr)

        C_r.append(corr.mean.real)
    print('C_r', jnp.array(C_r))

    ###### Fourier ######
    C_k = renormalized_corr(k, C_r)
    print('C_k', C_k)
    print('C_k_renorm', np.array(C_k)-C_r[0])
    print('numpy fft C_k', np.fft.fft(C_r))
    print('np.fft C_k_renorm', np.fft.fft(C_r)-C_r[0])

    ##############################################   
    ###### structure factor ######
    C_k = []
    for i in range(len(K)):
        C = corr_fourier(K[i], Nf, d,) # d neq rij !
        C =  vs.expect(C)
        C_k.append(C.mean.real)
    print('C_k tilde', np.array(C_k))
    print('C_k tilde_renorm', np.array(C_k)-C_r[0])

    ##############################################   
    ###### structure factor ######
    print()
    S_k = []
    for i in range(len(K)):

        s = structure_factor(K[i], Nf, rij,)
        s =  vs.expect(s)
        S_k.append(s.mean)
    print('structure factor', np.array(S_k))
    print('sf_renorm', np.array(S_k)-C_r[0])
    
    return energy_mean, standard, variance








R = g.positions
#print('R', R)
#print("Distance vectors:")
rij = distances_rij(R)
#print('distances', rij)

##################################################################################################################
def acceptance_callback(step_nr,log_data, driver):
    #if isinstance(driver.state, nk.vqs.ExactState) or not hasattr(driver.state.sampler_state, "acceptance"):
        #return True
    #else:
    acc = driver.state.sampler_state.acceptance
    if acc is None:
        acc = np.nan
    log_data["acceptance"] = float(acc)
    if mpi.node_number == 0:
        print("acc=", acc, flush=True)
    return True



###### symmetries ######

#### define symmetries: translations, rotations, reflexions,...
symm = g.translation_group()

character = symm.character_table()-0j
character = HashableArray(character[0])
print('character', jnp.array(character))


##################################################################################################################

R = g.positions

def distances_rij(r, box_size=None):
    N_sites = r.shape[0]
    distances = np.zeros((N_sites, N_sites, 2))

    for i in range(N_sites):
        for j in range(N_sites):
            diff = r[i] - r[j]

            if box_size is not None:
                # Apply periodic boundary conditions
                diff = diff - np.round(diff / box_size) * box_size

            distances[i, j] = diff
                        
    return distances


d = distances_rij(R, np.array([L,L]))


##################################################################################################################

####### optimize all parameters simultaneously #######
vs = run_ansatz(V)

####### validation: compute energy mean, error of mean and variance
energy = mean_en_var(ham, vs)
print('E', np.array(energy))
np.savetxt(f"output/en_mean_var_{V}.txt",[energy])
print('acceptance', vs.sampler_state.n_accepted_proc, vs.sampler_state.acceptance)
