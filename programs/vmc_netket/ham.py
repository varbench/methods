import netket as nk
from netket import experimental as nkx


def ColoredJ1J2(extent, pbc, *, back_diag=True):
    L1, L2 = extent

    def k(i, j):
        return (i % L1) * L2 + (j % L2)

    edges = []
    for i in range(L1):
        for j in range(L2 - 1 + pbc):
            edges.append((k(i, j), k(i, j + 1), 0))
    for i in range(L1 - 1 + pbc):
        for j in range(L2):
            edges.append((k(i + 1, j), k(i, j), 0))
    for i in range(L1 - 1 + pbc):
        for j in range(L2 - 1 + pbc):
            edges.append((k(i, j), k(i + 1, j + 1), 1))
            if back_diag:
                edges.append((k(i + 1, j), k(i, j + 1), 1))

    graph = nk.graph.Graph(edges=edges, n_nodes=L1 * L2)
    return graph


def Hubbard(hilbert, graph, U):
    def c(i, sz):
        return nkx.operator.fermion.destroy(hilbert, i, sz)

    def c_dag(i, sz):
        return nkx.operator.fermion.create(hilbert, i, sz)

    def n(i, sz):
        return nkx.operator.fermion.number(hilbert, i, sz)

    up = +1
    down = -1
    H = 0
    for i, j in graph.edges():
        for sz in [up, down]:
            H -= c_dag(i, sz) * c(j, sz) + c_dag(j, sz) * c(i, sz)
    for i in graph.nodes():
        H += U * n(i, up) * n(i, down)
    return H


def tVModel(hilbert, graph, V):
    def c(i):
        return nkx.operator.fermion.destroy(hilbert, i)

    def c_dag(i):
        return nkx.operator.fermion.create(hilbert, i)

    def n(i):
        return nkx.operator.fermion.number(hilbert, i)

    H = 0
    for i, j in graph.edges():
        H -= c_dag(i) * c(j) + c_dag(j) * c(i)
        H += V * n(i) * n(j)
    return H
