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


def to_ls(H: AbstractOperator, symmetries: ls.Symmetries) -> ls.Operator:
    if isinstance(H, SpecialHamiltonian):
        H = H.to_local_operator()

    if isinstance(H, LocalOperatorBase):
        return local_to_ls(H, symmetries)
    elif isinstance(H, FermionOperator2ndBase):
        return fermion_to_ls(H, symmetries)
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


def local_to_ls(H: LocalOperatorBase, symmetries: ls.Symmetries) -> ls.Operator:
    hilbert = H.hilbert
    assert isinstance(hilbert, nk.hilbert.Spin)
    if hilbert._total_sz is None:
        hamming_weight = None
    else:
        hamming_weight = hilbert.size // 2 + hilbert._total_sz * 2
    # TODO: Symmetries
    basis = ls.SpinBasis(
        hilbert.size,
        hamming_weight=hamming_weight,
        spin_inversion=1,
        symmetries=symmetries,
    )
    basis.build()

    exprs = None
    for matrix, acting_on in zip(H.operators, H.acting_on):
        paulis, weights = matrix_to_paulis(matrix)
        for pauli, weight in zip(paulis, weights):
            expr = pauli_to_expr(pauli, acting_on)
            if exprs is None:
                exprs = weight * ls.Expr(expr)
            else:
                exprs += weight * ls.Expr(expr)
    H = ls.Operator(basis, exprs)
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


def fermion_to_ls(H: FermionOperator2ndBase, symmetries: ls.Symmetries) -> ls.Operator:
    assert symmetries is None

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
        if exprs is None:
            exprs = weight * ls.Expr(expr)
        else:
            exprs += weight * ls.Expr(expr)
    H = ls.Operator(basis, exprs)
    return H


def get_symmetries(ham: str, extents: tuple[int, ...], pbc: bool) -> ls.Symmetries:
    ham_dim = len(extents)
    assert len(extents) == ham_dim

    symmetries = []
    if ham_dim == 1:
        N = extents[0]
        if ham in ["ising", "heis"]:
            if pbc:
                symmetries.append(ls.Symmetry([(i + 1) % N for i in range(N)], 0))
            symmetries.append(ls.Symmetry([N - i - 1 for i in range(N)], 0))
            symmetries = ls.Symmetries(symmetries)
            return symmetries
        else:
            return None
    elif ham_dim == 2:
        L1, L2 = extents
        N = L1 * L2

        def k(i, j):
            return (i % L1) * L2 + (j % L2)

        def ijrange():
            for k in range(N):
                yield divmod(k, L2)

        if ham in ["ising", "heis", "heis_tri", "j1j2"]:
            if pbc:
                symmetries.append(ls.Symmetry([k(i + 1, j) for i, j in ijrange()], 0))
                symmetries.append(ls.Symmetry([k(i, j + 1) for i, j in ijrange()], 0))
            if ham != "heis_tri":
                symmetries.append(
                    ls.Symmetry([k(L1 - i - 1, j) for i, j in ijrange()], 0)
                )
                symmetries.append(
                    ls.Symmetry([k(i, L2 - j - 1) for i, j in ijrange()], 0)
                )
            if L1 == L2:
                symmetries.append(ls.Symmetry([k(j, i) for i, j in ijrange()], 0))
            symmetries = ls.Symmetries(symmetries)
            return symmetries
        else:
            return None
    else:
        return None
