import numpy as np
import json
import os

def make_coulombintra(dir_def,file_coulombintra,Nsite,U):
    f = open(dir_def+"/"+file_coulombintra,"w")
    f.write("====\n")
    f.write("NCoulombIntra "+"{}\n".format(Nsite))
    f.write("====\n")
    f.write("==== CoulombIntra\n")
    f.write("====\n")
    for i in range(Nsite):
        f.write("{:5d} {:13.10f}\n".format(i,U))
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
    U = loaded_dict["U"]

    if not os.path.exists(dir_def):
        os.makedirs(dir_def)

    file_coulombintra = "coulombintra.def"
    make_coulombintra(dir_def,file_coulombintra,Nsite,U)
