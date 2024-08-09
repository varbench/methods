import numpy as np
import json
import os

def make_qptransidx(output,filename,Nsite,Lx,Ly):
    f = open(output+"/"+filename,"w")
    f.write("====\n")
    f.write("NQPTrans {}\n".format(Nsite))
    f.write("====\n")
    f.write("==== TrIdx_TrWeight_and_TrIdx_i_xi\n")
    f.write("====\n")
    for i in range(Nsite):
        f.write("{:5d} {:13.10f}\n".format(i,1.0))
    for i in range(Nsite):
        ix = i%Lx
        iy = i//Lx
        for j in range(Nsite):
            jx = j%Lx
            jy = j//Lx
            # P case
            kx = (ix+jx)%Lx
            ky = (iy+jy)%Ly
            k = kx + Lx*ky
            f.write("{:5d} {:5d} {:5d}\n".format(i,j,k))
            # AP case --> consider later
            #f.write("{:5d} {:5d} {:5d} {:5d}\n".format(i,j,k,1))
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

    if not os.path.exists(dir_def):
        os.makedirs(dir_def)

    file_qptransidx = "qptransidx.def"
    make_qptransidx(dir_def,file_qptransidx,Nsite,Lx,Ly)
