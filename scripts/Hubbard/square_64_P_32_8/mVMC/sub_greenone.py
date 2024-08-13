import numpy as np
import json
import os

def make_greenone(output,filename,Nsite):
    f = open(output+"/"+filename,"w")
    f.write("====\n")
    f.write("NCisAjs {}\n".format(2*Nsite))
    f.write("====\n")
    f.write("====\n")
    f.write("====\n")
    for i in range(Nsite):
        f.write("{:5d} 0 {:5d} 0\n".format(i,i))
        f.write("{:5d} 1 {:5d} 1\n".format(i,i))
    f.close()

if __name__ == "__main__":
    file_json = "dat_output_makedef.json"

    with open(file_json,"r") as json_file:
        loaded_dict = json.load(json_file)
    # print(loaded_dict)

    dir_def = loaded_dict["dir_def"]

    Lx = loaded_dict["Lx"]
    Ly = loaded_dict["Ly"]
    Nsite = Lx*Ly

    if not os.path.exists(dir_def):
        os.makedirs(dir_def)

    file_greenone = "greenone.def"
    make_greenone(dir_def,file_greenone,Nsite)
