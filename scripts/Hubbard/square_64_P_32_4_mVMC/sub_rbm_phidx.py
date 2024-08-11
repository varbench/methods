import numpy as np
import json
import os

def make_rbm(dir_def,file_rbm,Nsite,Lx,Ly,alpha):
    f = open(dir_def+"/"+file_rbm,"w")
    f.write("====\n")
    f.write("NRBM_HiddenLayerIdx "+"{}\n".format(2*alpha*Nsite))
    f.write("ComplexType 1\n")
    f.write("i s k RBM_PhysHidden_Idx\n")
    f.write("====\n")
    for j in range(alpha*Nsite):
        jalpha = j//Nsite
        jsite = j%Nsite
        jx = jsite%Lx
        jy = jsite//Lx
        for i in range(Nsite):
            ix = i%Lx
            iy = i//Lx
            for spin in range(2):
                kx = (ix+jx)%Lx
                ky = (iy+jy)%Lx
                k = 2*(kx+Lx*ky) + spin + 2*jalpha*Nsite
                f.write("{:d} {:d} {:d} {:d}\n".format(i,spin,j,k))
    for i in range(2*alpha*Nsite):
        f.write("{:d} {:d}\n".format(i,1))
    f.close()

if __name__ == "__main__":
    file_json = "dat_output_makedef.json"
    with open(file_json,"r") as json_file:
        loaded_dict = json.load(json_file)
    #print(loaded_dict)

    dir_def = loaded_dict["dir_def"]
    Lx = loaded_dict["Lx"]
    Ly = loaded_dict["Ly"]
    Nsite = Lx*Ly
    alpha = loaded_dict["alpha"]

    if not os.path.exists(dir_def):
        os.makedirs(dir_def)

    file_rbm = "rbm_phidx.def"
    make_rbm(dir_def,file_rbm,Nsite,Lx,Ly,alpha)
