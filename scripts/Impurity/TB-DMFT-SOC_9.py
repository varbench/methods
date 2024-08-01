import numpy as np
import scipy.linalg as lg
from scipy.optimize import minimize

from triqs.gf import *
import forktps as ftps
from forktps.solver import DMRGParams, TevoParams


from forktps.DiscreteBath import *
from forktps.Helpers import getX, MakeGFstruct
from h5 import *
from scipy.integrate import simps
from scipy.integrate import trapz


Nbath = 9
bondList = [310, 330, 350]

# interaction parameters
U = 2.3
J = 0.4
Up = U - 2.0 * J

withSOC = True
# read in Delta and e0
if withSOC:
    print("---- SoC calculations ........")
    with HDFArchive("SrRuO214_SOC_InFile.h5", "r") as ar:
        Delta_w = ar["Delta_w"].copy()
        e0_mat = ar["e0_mat"].copy()
else:
    print("---- NoSoC calculations ........")
    with HDFArchive("SrRuO214_NOSOC_InFile.h5", "r") as ar:
        Delta_w = ar["Delta_w"].copy()
        e0_mat = ar["e0_mat"].copy()

for bond in bondList:
    print("------ calculating Nb=", Nbath, ", bond=", bond)
    Hint = ftps.solver_core.HInt(u=U, j=J, up=Up, dd=False)
    S = ftps.Solver(gf_struct=MakeGFstruct(Delta_w), nw=3001, wmin=-10.0, wmax=10.0)

    # discritize the bath
    if withSOC:
        bath = DiscretizeBath(
            Delta=Delta_w,
            Nb=Nbath,
            PlaceAt={"ud_0": [0.0, 0.0, 0.0], "ud_1": [0.0, 0.0, 0.0]},
        )
    else:
        bath = DiscretizeBath(
            Delta=Delta_w,
            Nb=Nbath,
            PlaceAt={"up": [0.0, 0.0, 0.0], "down": [0.0, 0.0, 0.0]},
        )
    S.b = bath
    # fill e0
    e0 = ftps.solver_core.Hloc(MakeGFstruct(Delta_w), withSOC)
    for name, g in Delta_w:
        e0.Fill(name, e0_mat[name])
    S.e0 = e0
    ## solver params

    paramsDMRG = DMRGParams(
        maxmI=bond, maxmIB=bond, maxmB=bond, twI=1e-12, twIB=1e-12, twB=1e-12
    )

    # tevo params also contain other paramters as for example the time step dt and the number of time steps
    paramsTevo = TevoParams(
        maxmI=10,
        maxmIB=60,
        maxmB=100,
        twI=1e-9,
        twIB=1e-10,
        twB=1e-11,
        dt=0.1,
        time_steps=1,
    )

    S.solve(h_int=Hint, tevo=paramsTevo, params_GS=paramsDMRG, eta=0.1, calc_me=[])
    if withSOC:
        fileName = "Results_SrRuO214_SOC.h5"
    else:
        fileName = "Results_SrRuO214_NOSOC.h5"
    with HDFArchive(fileName, "a") as ar:
        ar["E0_Nbath_%s_bond_%s" % (Nbath, bond)] = S.E0
        ar["GSVariance_Nbath_%s_bond_%s" % (Nbath, bond)] = S.GSVariance
        ar["NPartGS_Nbath_%s_bond_%s" % (Nbath, bond)] = S.NPartGS
        ar["SectorEnergies_Nbath_%s_bond_%s" % (Nbath, bond)] = S.SectorEnergies
        ar["SectorImpOccs_Nbath_%s_bond_%s" % (Nbath, bond)] = S.SectorImpOccs
