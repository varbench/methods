import netket as nk
import flax.linen as nn
import jax
import jax.numpy as jnp
from jax.nn.initializers import lecun_normal, zeros
from netket.utils.types import NNInitFunc

from typing import Optional, Any, Union, Sequence, Tuple, Callable

from utils import _where_idx, _extract_cols, _log_det, _single_part

default_kernel_init = lecun_normal()
default_bias_init = zeros

class MeanBackflowSlater(nn.Module):

    L: int
    """The size of the system lattice is L^D. """
    D: int
    """The dimension of the system lattice. """
    Nf: int
    """The number of fermions on the discrete lattice. """
    Ns: int
    """The total number of sites on the lattice. """


    mf_orbitals: bool = True
    """Bool indicating whether to use mean field orbitals or plane waves. Defaults to mean field orbitals --> optimized matrix. """
    jastrow: bool = True
    """Bool indicating whether to use a (RBM) jastrow term which introduces correlations to our ansatz. """

    backflow: bool = True
    """Bool indicating whether to use a (MLP) backflow correction function. """

    param_dtype: Any = jnp.complex128
    """The dtype of the weights."""
    precision: Any = None
    """The numerical precision of the computation see :class:`jax.lax.Precision` for details."""

    ### MLP backflow params
    hidden_dims: Optional[Union[int, Tuple[int, ...]]] = None
    """The size of the hidden layers, excluding the output layer."""
    hidden_dims_alpha: Optional[Union[int, Tuple[int, ...]]] = None
    """The size of the hidden layers provided as number of times the input size.
    One must choose to either specify this or the hidden_dims keyword argument"""
    hidden_activations: Optional[Union[Callable, Tuple[Callable, ...]]] = nk.nn.reim_relu
    """The nonlinear activation function after each hidden layer.
    Can be provided as a single activation,
    where the same activation will be used for every layer."""
    output_activation: Optional[Callable] = None
    """The nonlinear activation at the output layer.
    If None is provided, the output layer will be essentially linear."""
    use_hidden_bias: bool = True
    """If True uses a bias in the hidden layer."""
    use_output_bias: bool = True
    """If True adds a bias to the output layer."""
    kernel_init: NNInitFunc = default_kernel_init
    """Initializer for the Dense layer matrix."""
    bias_init: NNInitFunc = default_bias_init
    """Initializer for the biases."""

    activation: Callable = nk.nn.reim_relu
    """The nonlinear activation function between layers."""

    #### RBM params
    activation_rbm: Callable = nk.nn.log_cosh
    """The nonlinear activation function between layers."""
    alpha: int = 1
    """The features density. Number of features equal to alpha * input.shape[-1]."""

    use_hidden_bias_rbm: bool = True
    """if True uses a bias in the dense layer (hidden layer bias)."""
    use_visible_bias_rbm: bool = True
    """if True adds a bias to the input not passed through the nonlinear layer."""

    kernel_init_rbm: NNInitFunc = default_kernel_init
    """Initializer for the Dense layer matrix."""
    hidden_bias_init_rbm: NNInitFunc = default_bias_init
    """Initializer for the hidden bias."""
    visible_bias_init_rbm: NNInitFunc = default_bias_init
    """Initializer for the visible bias."""




    # Backflow MLP architecture dim = = [# input nodes, # hidden layer nodes,...., # output nodes]

    #def setup(self):
        # MLP architecture - dim = = [# input nodes, # first hidden layer nodes,...., # output nodes]
        #self.hidden_layers = (self.alpha * self.Ns,)  # Ns is the total number of sites (= # inpute nodes)
        #self.out_dim = (self.Nf * self.Ns)

    @nn.compact
    def __call__(self, n):

        idx = jax.vmap(_where_idx, in_axes=(0, None))(n, self.Nf)[0]  # get idx of where particles lie given configuration |n>

        # Normalization of input: rescale the input to a signed binary notation of occupation number basis
        # centering the data around zero has shown to improve the NN performance
        n = (2*n-1)

        if self.mf_orbitals == True:
            # Mean field wave function --> optimize orbitals
            phi_j = self.param("lambda", default_kernel_init, (self.Nf, self.Ns), jnp.complex128) #jax.nn.initializers.normal(), (self.Nf, self.Ns), complex)  # var param matrix Nf * Nsites

        else:
            # Single particle orbitals
            phi_j = _single_part(self.L, self.D, self.Nf)


        if self.backflow == True:

            #bf = nk.models.MLP(output_dim=self.Nf*self.Ns,hidden_dims=2, param_dtype=self.param_dtype,
            #precision= self.precision, hidden_activations=self.activation)(n)

            #bf = jnp.expand_dims(bf, 0)
            #bf = bf.reshape(n.shape[0], self.Nf, self.Ns)
            bf = nn.Dense(2, param_dtype=self.param_dtype)(n)
            bf = jax.nn.tanh(bf)
            # last layer, outputs N x Nf values
            bf = nn.Dense(self.Nf*self.Ns, param_dtype=self.param_dtype)(bf)
            # reshape into M and add
            bf = jnp.expand_dims(bf, 0)
            bf = bf.reshape(n.shape[0], self.Nf, self.Ns)




            phi_bf = phi_j * bf
            phi = jax.vmap(_extract_cols, in_axes=(0, 0))(phi_bf,
                                                     idx)  # for each sample x take for each row in phi_param the Nf active indices (idx)

        else:
            phi = jax.vmap(_extract_cols, in_axes=(None, 0))(phi_j, idx) # for each sample x take for each row in phi_param the Nf active indices (idx)


        det = jax.vmap(_log_det)(phi)  # calculate (log) determinant of phi

        if self.jastrow == True:
            x = nk.models.RBM(alpha=self.alpha, precision=self.precision, use_hidden_bias=self.use_hidden_bias_rbm,
            use_visible_bias=self.use_visible_bias_rbm, param_dtype=self.param_dtype, activation = self.activation_rbm,
            kernel_init=self.kernel_init_rbm, hidden_bias_init=self.hidden_bias_init_rbm,
            visible_bias_init=self.visible_bias_init_rbm)(n)
            y = det + x
        else:
            y = det

        return y
