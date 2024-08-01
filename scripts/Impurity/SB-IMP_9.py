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


Nbath=9
U=2.0
bondSet=100
bondList = [bondSet]
Filling=1.0
model = "IMP"


# interaction parameters

# read in Delta and e0
print("---- one band calculations U = ", U,", N=", Filling, " ", model ,  "........")
with HDFArchive("OneBand_InFile.h5", 'r') as ar:
    Delta_w = ar["Delta_w_" + model ].copy()
    e0_mat = ar["e0_mat_" + model ].copy()
    print(e0_mat)


    
for bond in bondList:
    print("------ calculating Nb=", Nbath, ", bond=", bond)
    Hint = ftps.solver_core.HInt(u=U, j=0, up=0, dd=True)
    S = ftps.Solver(gf_struct=MakeGFstruct(Delta_w), nw=3001, wmin=-10.0, wmax=10.0)

    # discritize the bath
    if abs(U-8.0) < 0.1 and Filling==1.0:
        bath = DiscretizeBath(
            Delta=Delta_w,
            Nb=Nbath+1
            #gap=[-1.0,1.0] 
        )
    else:
        bath = DiscretizeBath(
            Delta=Delta_w,
            Nb=Nbath,
            PlaceAt={"up": [0.0], "dn": [0.0]},
        )
    S.b = bath
    # fill e0
    e0 = ftps.solver_core.Hloc(MakeGFstruct(Delta_w), False)
    for name, g in Delta_w:
        e0.Fill(name, e0_mat[name])
    S.e0 = e0
    tw = 1E-40 # truncated weight for all links
    DMRGSec = ftps.solver.DMRGParams( tw=1E-12, maxmI =100, maxmIB=100, maxmB=100, 
                                            nmax = 8, normErr = 1E-12, conv = 1E-14,
                                      sweeps =30, napph = 0, prep_napph=3)

    DMRG = ftps.solver.DMRGParams( tw=tw, maxm=bond,
                                   nmax = 16, normErr = 1E-14, conv = 1E-14,
                                   sweeps = 60, napph = 0, DMRGMethod="TwoSite" , prep_napph=5,
                                   prep_imagTevo=False, prep_dtau=0.02, prep_time_steps=200, prep_method="TEBD")

    tevo = ftps.solver.TevoParams(imag_tevo = False, dt = 0.01, time_steps = 1,
                                       nmax = 200, conv=1E-12, normErr=1E-12,
                                       tw=tw, maxmI =bond, maxmIB=300, maxmB=400)

    calc_me=["up", 0]
    S.solve(h_int=Hint, tevo=tevo,
            params_partSector = DMRGSec,
            params_GS = DMRG,
            eta=0.1, calc_me=[calc_me])
    
    fileName = "Results_OneBand.h5"
    with HDFArchive(fileName, 'a') as ar:
        ar["bath_Nbath_%s_bond_%s" % (Nbath, bond)] = bath
        ar["E0_Nbath_%s_bond_%s" % (Nbath, bond)] = S.E0
        ar["GSVariance_Nbath_%s_bond_%s" % (Nbath, bond)] = S.GSVariance
        ar["NPartGS_Nbath_%s_bond_%s" % (Nbath, bond)] = S.NPartGS
        ar["SectorEnergies_Nbath_%s_bond_%s" % (Nbath, bond)] = S.SectorEnergies
        ar["SectorImpOccs_Nbath_%s_bond_%s" % (Nbath, bond)] = S.SectorImpOccs
import os
os.remove("OneBand_InFile.h5")