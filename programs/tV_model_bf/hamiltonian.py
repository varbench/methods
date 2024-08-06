# import necessary libraries
import netket as nk
from netket import experimental as nkx


import jax.numpy as jnp
from scipy.sparse.linalg import eigsh


def t_v_model(L,D,t, U, Nf):

    # create the graph our fermions can hop on
    g = nk.graph.Hypercube(length=L, n_dim=D, pbc=True)
    n_sites = g.n_nodes

    # create a hilbert space
    hi = nkx.hilbert.SpinOrbitalFermions(n_sites,n_fermions=(Nf))

    # create an operator representing fermi hubbard interaction
    def cdag(site):
        return nkx.operator.fermion.create(hi, site)#nkx.operator.fermion.create(hi, site)


    def c(site):
        return nkx.operator.fermion.destroy(hi, site)


    def nc(site):
        return nkx.operator.fermion.number(hi, site)

    ham = 0.0

    for u, v in g.edges():
        ham += -t * cdag(u) * c(v) - t * cdag(v) * c(u)
    for u, v in g.edges():
        ham += U * nc(u) * nc(v)
    return hi, g, ham


def ED(ham):
    #exact = jnp.linalg.eigvalsh(ham.to_dense())
    exact = eigsh(ham.to_sparse(), k=2, which='SA')
    return exact
