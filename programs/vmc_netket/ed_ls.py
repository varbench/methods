#!/usr/bin/env python3
#
# SpinED (https://github.com/twesterhout/spin-ed) uses primme instead of
# scipy's eigsh but it's hard to build

from scipy.sparse.linalg import eigsh

from nk2ls import to_ls
from vmc import get_ham

_, _, H = get_ham()
H = to_ls(H)
energy = eigsh(H, k=1, which="SA", return_eigenvectors=False)
energy = energy[0]
print("energy", energy)
