import os
os.environ['OMP_NUM_THREADS'] = '4' # set number of OpenMP threads to run in parallel

import quspin

from quspin.operators import hamiltonian # operators
from quspin.basis import spinless_fermion_basis_general # spin basis constructor
import numpy as np # general math functions

from scipy.sparse.linalg import eigsh

#
###### define model parameters ######
Lx, Ly = 6, 6 # linear dimension of spin 1 2d lattice
Nf = 16# number of fermions
N_2d = Lx*Ly # number of sites for spin 1
#

import argparse

# Create the argument parser
parser = argparse.ArgumentParser()

# Add positional arguments
parser.add_argument("--V", type=float, help="V", required=True)

# Parse the arguments
args = parser.parse_args()

J=1. # hopping matrix element
U=args.V # onsite interaction
#mu=0.5 # chemical potential
#
###### setting up user-defined symmetry transformations for 2d lattice ######
s = np.arange(N_2d) # sites [0,1,2,....]
x = s%Lx # x positions for sites
y = s//Lx # y positions for sites
T_x = (x+1)%Lx + Lx*y # translation along x-direction
T_y = x +Lx*((y+1)%Ly) # translation along y-direction
P_x = x + Lx*(Ly-y-1) # reflection about x-axis
P_y = (Lx-x-1) + Lx*y # reflection about y-axis
#
###### setting up bases ######
basis_2d=spinless_fermion_basis_general(N_2d,Nf=Nf,kxblock=(T_x,0),kyblock=(T_y,1))#,pxblock=(P_x,0),pyblock=(P_y,0))
#
###### setting up hamiltonian ######
# setting up site-coupling lists
hopping_left=[[-J,i,T_x[i]] for i in range(N_2d)] + [[-J,i,T_y[i]] for i in range(N_2d)]
hopping_right=[[+J,i,T_x[i]] for i in range(N_2d)] + [[+J,i,T_y[i]] for i in range(N_2d)]
#potential=[[-mu,i] for i in range(N_2d)]
interaction=[[U,i,T_x[i]] for i in range(N_2d)] + [[U,i,T_y[i]] for i in range(N_2d)]
#
static=[["+-",hopping_left],["-+",hopping_right],["nn",interaction]]
# build hamiltonian
H=hamiltonian(static,[],basis=basis_2d,dtype=np.complex128)

print('H constructed!:)')
# diagonalise H
#E=H.eigvalsh()

print('now we calculate eigvals')
E,V = eigsh(H.aslinearoperator(),k=1,which="SA")


print(list(E))
print(len(list(E)))



