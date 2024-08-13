import numpy as np
import json
import os

def make_coulombinter(dir_def,file_coulombinter,Lx,Ly,Nsite,V,V2):
    f = open(dir_def+"/"+file_coulombinter,"w")
    f.write("====\n")
    if V2:
        n = 4
    else:
        n = 2
    f.write("NCoulombInter "+"{}\n".format(n*Nsite))
    #f.write("ComplexType 0\n")
    f.write("====\n")
    f.write("====\n")
    f.write("====\n")
    for i in range(Nsite):
        ix = i%Lx
        iy = i//Lx
        # hop right
        jx = (ix+1)%Lx
        jy = iy
        j = jx + Lx*jy
        f.write("{:5d} {:5d} {:13.10f}\n".format(i,j,V))
        # hop top
        jx = ix
        jy = (iy+1)%Ly
        j = jx + Lx*jy
        f.write("{:5d} {:5d} {:13.10f}\n".format(i,j,V))
    if V2:
        for i in range(Nsite):
            ix = i%Lx
            iy = i//Lx
            # hop right top
            jx = (ix+1)%Lx
            jy = (iy+1)%Ly
            j = jx + Lx*jy
            f.write("{:5d} {:5d} {:13.10f}\n".format(i,j,V2))
            # hop left top
            jx = (ix-1)%Lx
            jy = (iy+1)%Ly
            j = jx + Lx*jy
            f.write("{:5d} {:5d} {:13.10f}\n".format(i,j,V2))
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
    V = loaded_dict["V"]
    V2 = loaded_dict.get("V2", 0)

    if not os.path.exists(dir_def):
        os.makedirs(dir_def)

    file_coulombinter = "coulombinter.def"
    make_coulombinter(dir_def,file_coulombinter,Lx,Ly,Nsite,V,V2)
