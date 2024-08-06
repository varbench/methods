import numpy as np

import netket as nk
import jax
import jax.numpy as jnp
import flax.linen as nn

from jax.nn.initializers import lecun_normal, zeros
from flax.linen.dtypes import promote_dtype

from netket.jax import logsumexp_cplx
from jax.scipy.special import logsumexp
from netket.utils import HashableArray
from netket.utils.types import NNInitFunc, Array, DType
from netket.utils.group import PermutationGroup
from netket.graph import Lattice


from typing import Optional, Any, Union, Sequence, Tuple, Callable


from utils import _lattice_coord, _where_idx, _extract_cols, _log_det, _single_part, _extract, reorder_array
from symmetries import reshape_rotate, rotate, _symm_transf





default_kernel_init = lecun_normal()
default_bias_init = zeros


### Symmetrize: we want phi =  log sum_g exp(chi_g * det(gx) * symmRBM(gx))
# --> we need to define the determinants, a jastrow rbm and the corresponding symm characters

############ IMPORTANT: we re-define DenseSymm-layer which we use for the symmetric RBM, because we want to apply
#  the symm action g on the input x instead of the kernel, just to be consistent with how we
#  defined the symmetry transformation on the determinant: phi =  log sum_g exp(chi_g * det(gx) * symmRBM(gx))
# check in DenseSymmMatrix() below: x = jnp.take(x,jnp.asarray(self.symmetries), 2) and x = jnp.einsum('ijkl,mjl->imk',x,kernel)


##### Symmetric dense layer #####

class DenseSymmLayer(nn.Module):

    r"""Implements a symmetrized linear transformation over a permutation group
    using matrix multiplication."""

    symmetries: HashableArray
    """A group of symmetry operations (or array of permutation indices) over which the layer should be invariant.
        Numpy/Jax arrays must be wrapped into an :class:`netket.utils.HashableArray`.
    """
    features: int
    """The number of output features. Will be the second dimension of the output."""
    use_bias: bool = True
    """Whether to add a bias to the output (default: True)."""

    param_dtype: Any = jnp.complex128
    """The dtype of the weights."""
    kernel_init: NNInitFunc = default_kernel_init
    """Initializer for the kernel. Defaults to Lecun normal."""
    bias_init: NNInitFunc = default_bias_init
    """Initializer for the bias. Defaults to zero initialization."""


    def setup(self):
        self.symm = HashableArray(np.asarray(self.symmetries))
        self.n_symm, self.n_sites = np.asarray(self.symmetries).shape

    @nn.compact
    def __call__(self, x: Array) -> Array:
        """Applies the symmetrized linear transformation to the inputs along the last dimension.
        Args:
          x: The nd-array to be transformed. (needs to be 3dim)
        Returns:
          The transformed input.
        """

        # if x.ndim < 3, 'error'
        # TODO: Deprecated: Eventually remove and error if less than 3 dimensions
        if x.ndim < 3:
            old_shape = x.shape
            if x.ndim == 1:
                x = jnp.expand_dims(x, (0, 1))
            elif x.ndim == 2:
                x = jnp.expand_dims(x, 1)
            #symm_input_warning(old_shape, x.shape, "DenseSymm")

        in_features = x.shape[1]

        if self.use_bias:
            bias = self.param(
                "bias", self.bias_init, (self.features,), self.param_dtype)
        else:
            bias = None


        kernel = self.param("kernel",self.kernel_init,(self.features, in_features, self.n_sites),self.param_dtype,)


	    # apply symmetries to the input x instead of the kernel
        x = jnp.take(x,jnp.asarray(self.symm), 2)
        x = jnp.einsum('ijkl,mjl->imk',x,kernel) # outputs [batch, alpha, input features]


        if self.use_bias:
            # Convert symmetry-reduced bias of shape (features,) to the full bias of shape (..., features, 1).
            bias = jnp.expand_dims(bias, 1)


            x += bias



        return x


##### Symmetric RBM #####

class SymmRBM(nn.Module):

    r"""Constructs a symmetric RBM that consists of one dense symmetric layer and log_cosh activation function."""

    symmetries: Union[HashableArray, PermutationGroup]
    """A group of symmetry operations (or array of permutation indices) over which the layer should be invariant.
        Numpy/Jax arrays must be wrapped into an :class:`netket.utils.HashableArray`. """
    alpha: int = 1
    """The features density. Number of features equal to alpha * input.shape[-1]."""

    use_bias: bool = True
    """Whether to add a bias to the output (default: True)."""

    param_dtype: Any = jnp.complex128
    """The dtype of the weights."""

    kernel_init: NNInitFunc = default_kernel_init#lecun_normal()#default_equivariant_initializer
    """Initializer for the kernel. Defaults to Lecun normal."""
    bias_init: NNInitFunc = default_bias_init
    """Initializer for the bias. Defaults to zero initialization."""

    activation: Callable = nk.nn.log_cosh
    """The nonlinear activation function between layers."""

    def setup(self):
        self.n_symm, self.n_sites = np.asarray(self.symmetries).shape
        self.features = int(self.alpha * self.n_sites / self.n_symm)

    @nn.compact
    def __call__(self, n):

        # dense symmetric layer
        x = DenseSymmLayer(symmetries=self.symmetries, features=self.features*5, use_bias=self.use_bias, param_dtype=self.param_dtype,
        kernel_init=self.kernel_init, bias_init=default_bias_init)(n)

        # activation function
        x = self.activation(x)
        # summ along features axis and symm axis
        x = jnp.sum(x, axis=(-1,-2))
        return x


##### Symmetric dense neural network (g-CNN) #####
class GroupConvBackflow(nn.Module):
    out_dim: int
    hidden_layers: Tuple[int]

    Nf: int
    L: int
    symmetries: Union[HashableArray, PermutationGroup]
    activation: Callable = nk.nn.reim_relu
    last_bias: bool = True
    # last_linear: bool = True
    param_dtype: Any = complex
    alpha: int = 1


    def setup(self):

        self.n_symm, self.n_sites = np.asarray(self.symmetries).shape
        self.features = int(self.alpha * self.n_sites / self.n_symm)


    @nn.compact
    def __call__(self, n):

        x = DenseSymmLayer(symmetries=self.symmetries, features=self.features, use_bias=False, param_dtype=jnp.complex128,)(n)

        x = nk.nn.reim_relu(x)

        return x






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #










##### Convolutional neural network #####
# this neural network code was written by Jannes Nys as an alternative to flax.linen.Conv(features=.., kernel_size=(3, 3), strides=(1, 1), padding="CIRCULAR")(inputs)

class ConvFastCircular(nn.Module):
    """Fast Convolution Module for circular padding`.

    Attributes:
    features: number of convolution filters.
    kernel_size: shape of the convolutional kernel. For 1D convolution,
        the kernel size can be passed as an integer. For all other cases, it must
        be a sequence of integers.
    use_bias: whether to add a bias to the output (default: True).
    dtype: the dtype of the computation (default: infer from input and params).
    param_dtype: the dtype passed to parameter initializers (default: complex128).
    kernel_init: initializer for the convolutional kernel.
    bias_init: initializer for the bias.
    """
    features: int
    kernel_size: Sequence[int]
    use_bias: bool = True
    dtype: Optional[DType] = None
    param_dtype: DType = jnp.complex128
    kernel_init: NNInitFunc = default_kernel_init
    bias_init: NNInitFunc = zeros

    @nn.compact
    def __call__(self, inputs: Array) -> Array:
        """ Convolutional operation
        Can be made faster by considering individual operations (we are repeating quite a few probably)
        Expected input shape = [*, *extent, dim]
        Expected output shape = [*, *extent, features]
        """

        if isinstance(self.kernel_size, int):
            raise TypeError('Expected Conv kernel_size to be a'
                            ' tuple/list of integers (eg.: [3, 3]) but got'
                            f' {self.kernel_size}.')
        else:
            kernel_size = tuple(self.kernel_size)

        # Combine all input batch dimensions into a single leading batch axis.
        num_batch_dimensions = inputs.ndim - (len(kernel_size) + 1)
        input_batch_shape = inputs.shape[:num_batch_dimensions]
        input_features = inputs.shape[-1]
        grid_shape = inputs.shape[num_batch_dimensions:-1]



        kernel_size_arr = np.array(kernel_size)

        #assert np.all(kernel_size_arr // 2 == 1), "kernels must be odd"
        max_shifts = (kernel_size_arr-1) // 2
        shifts = np.meshgrid(*[np.arange(-m, m+1) for m in max_shifts])
        shifts = list(map(np.ravel, shifts))

        n_shifts = shifts[0].shape[0]

        axes = tuple(np.arange(num_batch_dimensions, inputs.ndim - 1))
        shifted_input = [jnp.roll(inputs, sh, axis=axes) for sh in zip(*shifts)] # circular padding

        shifted_input = jnp.stack(shifted_input, axis=-2)

        shifted_input = shifted_input.reshape(*input_batch_shape, *grid_shape, n_shifts*input_features)

        # change last dimension only
        out = nn.Dense(self.features, use_bias=self.use_bias, dtype=self.dtype, param_dtype=self.param_dtype,
        kernel_init=self.kernel_init, bias_init=self.bias_init)(shifted_input)

        return out


##### CNN Backflow #####

def reshape_rotate(x, L):
    x = x.reshape(L,L)
    x = jnp.rot90(x, k=2)
    return jnp.flip(x, 0)


def rotate(x, L):

    #x = jnp.flip(x, 0)
    x = jnp.rot90(x, k=2)
    return x

class ConvBackflow(nn.Module):

    L: int
    """The size of the system lattice is L^D. """
    Ns: int
    """The total number of sites on the lattice. """

    Nf: int
    """The number of fermions on the discrete lattice. """

    symmetries: Union[HashableArray, PermutationGroup]
    """A group of symmetry operations (or array of permutation indices) over which the layer should be invariant.
        Numpy/Jax arrays must be wrapped into an :class:`netket.utils.HashableArray`. """
    activation: Callable = nk.nn.reim_relu
    """The nonlinear activation function between layers."""
    depth: int=1
    """The depth of the CNN backflow. """
    features: int = 1
    """The number of channels of the CNN backflow."""
    dtype: Optional[DType] = None
    """The dtype of the computation."""
    param_dtype: Any = jnp.complex128
    """The dtype of the weights."""

    kernel_size: Sequence[int] = (3,3)
    """The size of the CNN kernel."""

    use_bias: bool = True
    """Whether to add a bias to the output (default: True)."""

    kernel_init: NNInitFunc = default_kernel_init
    """Initializer for the kernel. Defaults to Lecun normal."""
    bias_init: NNInitFunc = default_bias_init
    """Initializer for the bias. Defaults to zero initialization."""

    def setup(self):
        layers = []
        for i in range(self.depth):
            layers.append(ConvFastCircular(features=self.features, kernel_size=self.kernel_size, use_bias=self.use_bias, dtype=self.dtype,
            param_dtype=self.param_dtype, kernel_init=self.kernel_init, bias_init=self.bias_init))
            #layers.append(nn.Conv(features=self.features, kernel_size=(3, 3), strides=(1, 1), padding="CIRCULAR", param_dtype = jnp.complex128))
        self.layers = layers
        self.layernorm = nn.LayerNorm(
            dtype=jnp.complex128, use_bias=True, use_scale=False, param_dtype = self.param_dtype
        )


    @nn.compact
    def __call__(self, n):

        n = n.reshape(n.shape[0], self.L,self.L,1)
        x = n

        residual = x.copy()

        if self.depth == 1:
            for i, layer in enumerate(self.layers):
                x = layer(x)
                x = self.activation(x)

        else:
            for i, layer in enumerate(self.layers):

                if i:
                    x = self.layernorm(x)
                    x = self.activation(x)

                x = layer(x)
                if i % 2:
                    x += residual
                    residual = x.copy()


        x = ConvFastCircular(self.Ns, (1,1),param_dtype = jnp.complex128)(x)
        #x  =nn.Conv(features=1, kernel_size=(1, 1), strides=(1, 1), padding="CIRCULAR",param_dtype = jnp.complex128)(x)

        x = x.reshape(x.shape[0],self.L,self.L, self.Ns)

        x = x.reshape(x.shape[0],self.L*self.L, self.Ns)

        t = np.asarray(self.symmetries)[:,0]
        x = jax.vmap(reorder_array, in_axes=(0,None))(x,t)

        return x



##### Wavefunction ansatz: symmetric Slater Backflow Jastrow #####
##### Symmetric dense neural network (g-CNN) #####
class DenseSymmBackflow(nn.Module):
    out_dim: int
    symmetries: Union[HashableArray, PermutationGroup]
    activation: Callable = nk.nn.reim_relu
    param_dtype: Any = complex
    alpha: int = 1

    kernel_init: NNInitFunc = default_kernel_init
    """Initializer for the kernel. Defaults to Lecun normal."""
    bias_init: NNInitFunc = zeros
    """Initializer for the bias. Defaults to zero initialization."""



    def setup(self):

        self.n_symm, self.n_sites = np.asarray(self.symmetries).shape
        self.features = int(self.alpha * self.n_sites / self.n_symm)


    @nn.compact
    def __call__(self, n):

        x = nk.nn.DenseSymm(symmetries=self.symmetries, features=self.features, use_bias=True,
                            param_dtype=self.param_dtype,kernel_init=self.kernel_init)(n)

        x = self.activation(x)

        x = nk.nn.DenseEquivariant(symmetries=self.symmetries, mode = "auto",features=self.features,
                                   use_bias=True, param_dtype=jnp.complex128,)(x)
        x = nk.nn.DenseEquivariant(symmetries=self.symmetries, mode = "auto",features=self.out_dim,
                                   use_bias=True, param_dtype=jnp.complex128,)(x)

        x = logsumexp_cplx(x, axis=(-1))
        return x


class JastrowTranslationInv(nn.Module):
    "Two-body Jastrow factor that is invariant under translations."

    r'''
    \displaystyle\log \psi(z) = 1/2 \sum_{ij} z_i W_{d(ij)} z_j
    '''
    graph : Lattice
    """The graph of the model."""

    n_neigh : int
    """The number of possible neighbors for site i on the lattice. For a periodic square lattice it is constant."""

    param_dtype : Any = jnp.complex128
    """The dtype of the weights."""

    jastrow_init : NNInitFunc = nn.initializers.normal() #nn.initializers.constant(np.array([1,2,3,4]))
    """Initializer for the jastrow parameters."""

    @nn.compact
    def __call__(self, n):

        # Jastrow variational parameters to be optimized
        J = self.param('J', self.jastrow_init, (self.n_neigh,), self.param_dtype)


        # Nearest-neighbor correlations
        corr = 0.5*jnp.einsum( '...i,ij,...j',n,J[self.graph.distances()],n )

        return corr


##### Wavefunction ansatz: symmetric Slater Backflow Jastrow #####

class SymmMeanBackflowSlater(nn.Module):

    L: int
    """The size of the system lattice is L^D. """
    D: int
    """The dimension of the system lattice. """
    Nf: int
    """The number of fermions on the discrete lattice. """
    Ns: int
    """The total number of sites on the lattice. """

    symmetries: Union[HashableArray, PermutationGroup]
    """A group of symmetry operations (or array of permutation indices) over which the layer should be invariant.
        Numpy/Jax arrays must be wrapped into an :class:`netket.utils.HashableArray`. """

    character: HashableArray
    """The character table corresponding to the given symmetry group. """

    graph : Lattice
    """The graph of the model."""

    mf_orbitals: bool = True
    """Bool indicating whether to use mean field orbitals or plane waves. Defaults to mean field orbitals --> optimized matrix. """
    jastrow: bool = True
    """Bool indicating whether to use a jastrow term which introduces correlations to our ansatz. """
    jastrow_rbm: bool = False
    """Bool indicating whether to use a SymmRBM jastrow term. If False, we use a two-body translation invariant Jastrow term. """



    backflow: bool = True
    """Bool indicating whether to use a (CNN) backflow correction function. """

    #### backflow params
    depth: int=1
    """The depth of the CNN backflow. """
    features: int = 1
    """The number of features of the NN backflow. """
    activation: Callable = nk.nn.reim_relu
    """The nonlinear activation function between layers."""
    kernel_size: Sequence[int] = (3,3)
    """The size of the CNN kernel."""
    use_bias: bool = True
    """Whether to add a bias to the output (default: True)."""
    kernel_init: NNInitFunc = default_kernel_init
    """Initializer for the kernel. Defaults to Lecun normal."""
    bias_init: NNInitFunc = default_bias_init
    """Initializer for the bias. Defaults to zero initialization."""

    #### RBM params

    alpha: int = 1
    """The features density. Number of features equal to alpha * input.shape[-1]."""
    use_bias_rbm: bool = True
    """Whether to add a bias to the output (default: True)."""
    kernel_init_rbm: NNInitFunc = default_kernel_init#lecun_normal()#default_equivariant_initializer
    """Initializer for the kernel. Defaults to Lecun normal."""
    bias_init_rbm: NNInitFunc = default_bias_init
    """Initializer for the bias. Defaults to zero initialization."""
    activation_rbm: Callable = nk.nn.log_cosh
    """The nonlinear activation function between layers."""



    dtype: Optional[DType] = None
    """The dtype of the computation."""
    param_dtype: Any = jnp.complex128
    """The dtype of the weights."""

    # Ns is the total number of sites (= # inpute nodes)
    last_bias: bool = True

    def setup(self):
        self.n_symm, self.n_sites = np.asarray(self.symmetries).shape




    @nn.compact
    def __call__(self, n):


        # get indices of where particles lie given configuration n
        idx = jax.vmap(_where_idx, in_axes=(0, None))(n, self.Nf)

        # apply symmetry transform into the indices
        symm_inv = np.asarray(self.symmetries)[self.symmetries.inverse] # we want g^{-1}^{-1} = g
        idx_g = jax.vmap(jax.vmap(_extract,in_axes=(0, None)), in_axes=(None, 0))(symm_inv, jnp.asarray(idx[0]))


        # Normalization of input: rescale the input to a signed binary notation of occupation number basis
        # centering the data around zero has shown to improve the NN performance
        n = (2*n-1)

        if self.mf_orbitals == True:
            # Mean field wave function --> optimize orbitals
            phi_j = self.param("phi_mf", default_kernel_init, (self.Nf, self.Ns), jnp.complex128)  # var param matrix of size Nf * Nsites


        else:
            # Single particle orbitals
            phi_j = _single_part(self.L, self.D, self.Nf)
            #print('phi_j plane waves', phi_j)


        # CNN backflow
        if self.backflow == True:
            # construct CNN backflow, we want bf to be of dimension number of orbitals (in our case = number of fermions) N_orbitals * number of sites N_sites
            #CNNBackflow vs CNNBackflow
            ################################################################################################################################
            # create a list that contains N_orbitals times the CNN backflow model (function)
            #bf_mu = [DenseSymmBackflow(out_dim=self.Ns,  symmetries=self.symmetries, activation=self.activation, param_dtype=self.param_dtype,
            #kernel_init=self.kernel_init, bias_init=self.bias_init) for i in range(self.Nf)]


            b = nn.vmap(ConvBackflow,in_axes=None, axis_size=self.Nf,    variable_axes={'params': 0},split_rngs={'params': True})

            bf_out = b(L=self.L, Ns = self.Ns, Nf=self.Nf, symmetries=self.symmetries, activation=self.activation, depth=self.depth,
            features=self.features, dtype=self.dtype, param_dtype=self.param_dtype, kernel_size=self.kernel_size, use_bias=self.use_bias,
            kernel_init=self.kernel_init, bias_init=self.bias_init)(n)
            bf_out = bf_out.transpose(1,2,0,3)

            #################################################################################################################################

            # for each sample n take for each row in phi_j the Nf active indices (idx)
            bf = jax.vmap(jax.vmap(_extract_cols, in_axes=(0, 0)), in_axes=(0,0))(bf_out,jnp.asarray(idx_g))

            #print('phi extracted', phi)

            # for each sample n take for each row in phi_j the Nf active indices (idx)
            phi_ = jax.vmap(_extract_cols, in_axes=(None, 0))(phi_j,jnp.asarray(idx_g))
            phi_ = phi_.transpose(0,2,1,3)
            # multiply mean field orbitals with backflow orbitals
            phi = phi_ * (1+bf)



        else:
            # symmetrize phi_j wrt to self.symmetries
            # for each sample n take for each row in phi_j the Nf active indices (idx)
            phi = jax.vmap(_extract_cols, in_axes=(None, 0))(phi_j,jnp.asarray(idx_g))
            phi = phi.transpose(0,2,1,3)

        # calculate (log) determinant of phi
        det = jax.vmap(_log_det)(phi.reshape(phi.shape[0], phi.shape[1],self.Nf, self.Nf))

        # jastrow correlation function
        if self.jastrow == True:
            if self.jastrow_rbm == True:
                # RBM jastrow consists of a dense symmetric layer
                x = SymmRBM(symmetries=self.symmetries, alpha=self.alpha, use_bias=self.use_bias_rbm, param_dtype=self.param_dtype,
                kernel_init = self.kernel_init_rbm, bias_init = self.bias_init_rbm, activation=self.activation_rbm)(n)
            else:
                x = JastrowTranslationInv(graph=self.graph,n_neigh=self.L)(n)

            def add(a,b):
                return a+b
            # add log slater determinant with log jastrow
            y = jax.vmap(add)(det,x)
        else:
            # no jastrow
            y = det

        # since we symmetrized our ansatz we need to introduce the characters associated to each symmetry transformation
        # symmetrize the wavefunction by averaging over all the symmetry transformations with logsumexp
        res = logsumexp_cplx(jnp.array(y),axis=(1), b=jnp.asarray(self.character).conj())
        return res
