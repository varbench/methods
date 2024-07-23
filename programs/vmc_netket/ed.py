#!/usr/bin/env python3

import netket as nk

from vmc import get_ham

_, _, H = get_ham()
energy = nk.exact.lanczos_ed(H, k=1)
energy = energy[0]
print("energy", energy)
