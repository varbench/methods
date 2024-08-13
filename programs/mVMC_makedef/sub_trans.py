import numpy as np
import json
import os

def make_trans(dir_def,file_trans,Lx,Ly,Nsite,t2):
    f = open(dir_def+"/"+file_trans,"w")
    f.write("====\n")
    if t2:
        n = 16
    else:
        n = 8
    f.write("NTransfer "+"{}\n".format(n*Nsite))
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
        f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(i,0,j,0,1.0,0))
        f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(j,0,i,0,1.0,0))
        f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(i,1,j,1,1.0,0))
        f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(j,1,i,1,1.0,0))
        # hop top
        jx = ix
        jy = (iy+1)%Ly
        j = jx + Lx*jy
        f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(i,0,j,0,1.0,0))
        f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(j,0,i,0,1.0,0))
        f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(i,1,j,1,1.0,0))
        f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(j,1,i,1,1.0,0))
    if t2:
        for i in range(Nsite):
            ix = i%Lx
            iy = i//Lx
            # hop right top
            jx = (ix+1)%Lx
            jy = (iy+1)%Ly
            j = jx + Lx*jy
            f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(i,0,j,0,t2,0))
            f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(j,0,i,0,t2,0))
            f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(i,1,j,1,t2,0))
            f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(j,1,i,1,t2,0))
            # hop left top
            jx = (ix-1)%Lx
            jy = (iy+1)%Ly
            j = jx + Lx*jy
            f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(i,0,j,0,t2,0))
            f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(j,0,i,0,t2,0))
            f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(i,1,j,1,t2,0))
            f.write("{:5d} {:5d} {:5d} {:5d} {:13.10f} {:.1f}\n".format(j,1,i,1,t2,0))
    # consider AP case later
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
    flag_ap = loaded_dict["APFlag"]
    t2 = loaded_dict.get("t2", 0)

    if not os.path.exists(dir_def):
        os.makedirs(dir_def)

    file_trans = "trans.def"
    make_trans(dir_def,file_trans,Lx,Ly,Nsite,t2)
