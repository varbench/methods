import jax
import numpy as np
import jax.numpy as jnp


def k_orbitals(L,Nf):
    k = 0
    kx = jnp.array([((2*i)*jnp.pi/L,(-(2*i)*jnp.pi/L)) for i in range(1,(L//2)+1)]).ravel() # coord in k space in x (y) direction
    kx = jnp.insert(kx,0,0)
    kx = kx[0:L]

    k = jnp.array([(i,j) for i in kx for j in kx]) # coordinates of the orbitals in k-space
    k = sort_k(k)
    k = k[0:Nf] # only take Nf orbitals
    return k

def sort_k(k):
    idx = jnp.argsort(jnp.linalg.norm(k,axis=1))
    k_sorted = k[idx]
    return k_sorted


def k_vector(L):
    kx = jnp.array([((2*i)*jnp.pi/L,(-(2*i)*jnp.pi/L)) for i in range(1,(L//2)+1)]).ravel() # coord in k space in x (y) direction
    kx = jnp.insert(kx,0,0)
    k = jnp.array([(i,j) for i in kx for j in kx]) # coordinates of the orbitals in k-space
    k = sort_k(k)
    
    return k
    
### PBC
def _lattice_coord(L,D,Nf): # defines the lattice coordinates of the sites in real and k-space
    r = jnp.array([(i,j) for i in range(L) for j in range(L)]) # coordinates of the lattice sites in real space
    k = k_orbitals(L,Nf) # coordinates of the lattice orbitals in k space
    return r, k




def _single_part(L,D,Nf): # return single part orbitals 

    coord = _lattice_coord(L,D,Nf)
    r = coord[0] # coordinates of the lattice sites in real space
    k = coord[1]
    phi = jnp.exp(1j*jnp.einsum('ij,kj->ik',k,r)) # take exp(ikr) as single spin orbital basis functions
    return phi



def _where_idx(x,Nf): # get idx of where particles lie given configuration |x>
    return jnp.where(x,size=Nf)

def _extract_cols(phi,idx):
    phi = phi[:,idx]
    return phi


def _log_det(phi): # return log slater determinant
    det = jnp.linalg.det(phi)
    log_det = jnp.linalg.slogdet(phi)
    return jnp.asarray(jnp.log(log_det[0])+log_det[1])


def get_max_conn(op, verbose=True, change=True, safety_factor=1.5):#safety_factor=1.1
    hi = op.hilbert
    _, mels = op.get_conn_padded(hi.random_state(key=jax.random.PRNGKey(0), size=(64*1024)))
    mels = mels.reshape(-1, mels.shape[-1])
    
    n_conn = (1-np.isclose(mels, 0)).sum(axis=-1)
    
    max_conn = np.max(n_conn)
    
    if change:
        assert safety_factor >= 1
        op._max_conn_size = int(max_conn*safety_factor)
        op._max_conn_size
        
    return op._max_conn_size


def distances_rij(r):
    N_sites = r.shape[0]
    distances = np.zeros((N_sites,N_sites, 2))

    for i in range(N_sites):
        for j in range(N_sites):
            distances[i, j] = r[i] - r[j]
            
    return distances


def reorder_array(batch, t):
    reordered_batch = jnp.take(batch, t, axis=0)
    return reordered_batch

def _extract(a,idx):
    a = a[idx]
    return a

