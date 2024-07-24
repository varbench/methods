import lattice_symmetries as ls
import netket as nk


# '0' -> '₀', '1' -> '₁', ...
def subscript_number(n):
    return "".join([chr(ord(c) + 8272) for c in str(n)])


def to_ls(H):
    if isinstance(H, nk.experimental.operator.FermionOperator2nd):
        return fermion_to_ls(H)
    else:
        raise TypeError(f"Unsupported operator {type(H)}: {H}")


def fermion_term_to_expr(term, spinful, N):
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


def fermion_to_ls(H):
    assert isinstance(H, nk.experimental.operator.FermionOperator2nd)

    hilbert = H.hilbert
    if hilbert.n_spin_subsectors == 1:
        basis = ls.SpinlessFermionBasis(hilbert.n_orbitals, hilbert.n_fermions)
        spinful = False
    else:
        assert hilbert.n_spin_subsectors == 2
        assert hilbert.n_fermions_per_spin[1] == hilbert.n_fermions_per_spin[0]
        # TODO: Set number_particles = Tuple[int, int] is documented but fails
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

    # import numpy as np
    # for i in range(hilbert.size):
    #     x = np.zeros(hilbert.size)
    #     x[i] = 1
    #     x = H @ x
    #     x[abs(x) < 1e-7] = 0
    #     print(x)

    return H
