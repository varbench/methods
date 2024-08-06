import jax
import jax.numpy as jnp


# translate r position in x, y and xy direction
def _translate_r(a):
    y = jnp.vstack([jnp.roll(jnp.roll(a, -j, axis=0), -i, axis=1).flatten() for i in range(a.shape[0])]
                           for j in range(a.shape[1]))
    return y


# symmetrize phi
def _symm_transf(phi, symm):
    phi_symm = jnp.take(phi,jnp.asarray(symm),1)
    phi_symm = jnp.asarray(phi_symm)
    return phi_symm


def _character(K,g): # character computed wrt momentum, K is the total momentum and g is the symmetry transformation vector, e.g. g = T_x, translation in x 
    return jnp.exp(-1j * jnp.sum(g*K, axis=1))


# identify active indices (particle positions) given configuration, e.g for config n = [0,1,0,1] indices 1 and 3 are active
def _extract(a,idx):
    a = a[idx]
    return a

def _parity(p):
    return jnp.linalg.det(jax.jacobian(jnp.sort)(p.astype(float))).astype(int)

# compute parities
def _parity_factor(symm, idx):

    # symmetry operations (or array of permutation indices)
    symm = jnp.asarray(symm)
    # active indices (indices of particle positions)
    idx = jnp.asarray(idx)
    
    # permute the active indices with the symmetry operations symm
    symm_extra = jax.vmap(jax.vmap(_extract, in_axes=(0,None)), in_axes=(None, 0))(symm, idx)  

    # parity factors
    parities =  jax.vmap(jax.vmap(_parity))(symm_extra) 

    return parities


# multiply a by b
def _mult(a, b):
    return a*b

def _add(a, b):
    return a+b


def reshape_rotate(x, L):
    x = x.reshape(L,L)
    return jnp.rot90(x, k=1)

def rotate(x, L):
    return jnp.rot90(x, k=2)#jnp.transpose(x)#k=2 works!

