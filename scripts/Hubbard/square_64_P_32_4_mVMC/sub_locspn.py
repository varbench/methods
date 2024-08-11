import numpy as np
import json
import os

def make_locspn(output,filename,Nsite):
    f = open(output+"/"+filename,"w")
    f.write("====\n")
    f.write("NlocalSpin 0\n")
    f.write("====\n")
    f.write("====\n")
    f.write("====\n")
    for i in range(Nsite):
        f.write("{:5d} {:5d}\n".format(i,0))
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

    if not os.path.exists(dir_def):
        os.makedirs(dir_def)

    file_locspn = "locspn.def"
    make_locspn(dir_def,file_locspn,Nsite)
