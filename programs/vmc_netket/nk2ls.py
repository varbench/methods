import itertools
from typing import Iterable

import lattice_symmetries as ls
import netket as nk
import numpy as np
from netket.experimental.operator._fermion_operator_2nd_base import (
    FermionOperator2ndBase,
)
from netket.operator import AbstractOperator
from netket.operator._hamiltonian import SpecialHamiltonian
from netket.operator._local_operator import LocalOperatorBase


# '0' -> '₀', '1' -> '₁', ...
def subscript_number(n: int) -> str:
    return "".join([chr(ord(c) + 8272) for c in str(n)])


def to_ls(H: AbstractOperator) -> ls.Operator:
    if isinstance(H, SpecialHamiltonian):
        H = H.to_local_operator()

    if isinstance(H, LocalOperatorBase):
        return local_to_ls(H)
    elif isinstance(H, FermionOperator2ndBase):
        return fermion_to_ls(H)
    else:
        raise TypeError(f"Unsupported operator {type(H)}: {H}")


def kron_prod(xs: Iterable[np.ndarray]):
    out = 1
    for x in xs:
        out = np.kron(out, x)
    return out


def matrix_to_paulis(
    matrix: np.ndarray, *, cutoff: float = 1e-10
) -> tuple[list[str], list[float | complex]]:
    I = np.array([[1, 0], [0, 1]], dtype=np.float64)
    X = np.array([[0, 1], [1, 0]], dtype=np.float64)
    Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    Z = np.array([[1, 0], [0, -1]], dtype=np.float64)

    assert matrix.ndim == 2
    assert matrix.shape[0] == matrix.shape[1]
    # Number of Pauli operators
    N = int(np.log2(matrix.shape[0]))
    assert 2**N == matrix.shape[0]
    # Can only handle small matrix
    assert N <= 10

    paulis_all = ["".join(x) for x in itertools.product("IXYZ", repeat=N)]
    vectors_all = [
        kron_prod(x).flatten() for x in itertools.product([I, X, Y, Z], repeat=N)
    ]
    target_vector = matrix.flatten()

    paulis = []
    weights = []
    for pauli, vector in zip(paulis_all, vectors_all):
        weight = (vector.conj() @ target_vector).item() / 2**N
        if abs(weight) < cutoff:
            continue
        if abs(weight.imag) < cutoff:
            weight = weight.real
        paulis.append(pauli)
        weights.append(weight)
    return paulis, weights


# TODO: Support constant (all I)
def pauli_to_expr(pauli: str, acting_on: tuple[int, ...]) -> str:
    assert len(pauli) == len(acting_on)
    expr = ""
    for c, i in zip(pauli, acting_on):
        if c == "I":
            continue

        if expr:
            expr += " "

        if c == "X":
            expr += "σˣ"
        elif c == "Y":
            expr += "σʸ"
        elif c == "Z":
            expr += "σᶻ"
        else:
            raise ValueError(f"Unsupported operator: {c}")

        expr += subscript_number(i)
    return expr


def local_to_ls(H: LocalOperatorBase) -> ls.Operator:
    hilbert = H.hilbert
    assert isinstance(hilbert, nk.hilbert.Spin)
    if hilbert._total_sz is None:
        hamming_weight = None
    else:
        hamming_weight = hilbert.size // 2 + hilbert._total_sz * 2
    # TODO: Symmetries
    basis = ls.SpinBasis(hilbert.size, hamming_weight=hamming_weight, spin_inversion=1)
    basis.build()

    exprs = None
    for matrix, acting_on in zip(H.operators, H.acting_on):
        paulis, weights = matrix_to_paulis(matrix)
        for pauli, weight in zip(paulis, weights):
            expr = pauli_to_expr(pauli, acting_on)
            # print(weight, expr)
            if exprs is None:
                exprs = weight * ls.Expr(expr)
            else:
                exprs += weight * ls.Expr(expr)
    H = ls.Operator(basis, exprs)

    # for i in range(2**hilbert.size):
    #     x = np.zeros(2**hilbert.size)
    #     x[i] = 1
    #     x = H @ x
    #     x[abs(x) < 1e-7] = 0
    #     print(x)

    return H


def fermion_term_to_expr(term: list[tuple[int, int]], spinful: bool, N: int) -> str:
    expr = ""
    for i, x in term:
        if expr:
            expr += " "

        if x == 0:
            expr += "c"
        else:
            expr += "c†"

        # TODO: The documented format is c†₀↑ but the parser actually accepts c†↑₀
        if spinful:
            if i < N:
                expr += "↓"
            else:
                i -= N
                expr += "↑"

        expr += subscript_number(i)
    return expr


def fermion_to_ls(H: FermionOperator2ndBase) -> ls.Operator:
    hilbert = H.hilbert
    if hilbert.n_spin_subsectors == 1:
        basis = ls.SpinlessFermionBasis(hilbert.n_orbitals, hilbert.n_fermions)
        spinful = False
    else:
        assert hilbert.n_spin_subsectors == 2
        assert hilbert.n_fermions_per_spin[1] == hilbert.n_fermions_per_spin[0]
        # TODO: Set number_particles = tuple[int, int] is documented but fails
        basis = ls.SpinfulFermionBasis(hilbert.n_orbitals, hilbert.n_fermions)
        spinful = True
    basis.build()

    exprs = None
    for term, weight in zip(H.terms, H.weights):
        expr = fermion_term_to_expr(term, spinful, hilbert.n_orbitals)
        # print(weight, expr)
        if exprs is None:
            exprs = weight * ls.Expr(expr)
        else:
            exprs += weight * ls.Expr(expr)
    H = ls.Operator(basis, exprs)
    return H
