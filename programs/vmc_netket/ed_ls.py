#!/usr/bin/env python3
#
# SpinED (https://github.com/twesterhout/spin-ed) uses primme instead of
# scipy's eigsh but it's hard to build

from scipy.sparse.linalg import eigsh

from args import args
from nk2ls import get_symmetries, to_ls
from vmc import get_ham

_, _, H = get_ham()

if args.ham_dim == 1:
    extents = (args.L,)
else:
    extents = (args.L, args.L2)

if args.boundary == "peri":
    pbc = True
elif args.boundary == "open":
    pbc = False
else:
    raise ValueError(f"Unknown boundary: {args.boundary}")

symmetries = get_symmetries(args.ham, extents, pbc)

H = to_ls(H, symmetries)
energy = eigsh(H, k=1, which="SA", return_eigenvectors=False)
energy = energy[0]
print("energy", energy)
