import numpy as np
import json
import os

def make_orbitalidx(output,filename,Nsite):
    f = open(output+"/"+filename,"w")
    f.write("====\n")
    f.write("NOrbitalIdx {}\n".format(Nsite**2))
    #f.write("ComplexType 0\n")
    f.write("ComplexType 1\n")
    f.write("====\n")
    f.write("====\n")
    cnt = 0
    for i in range(Nsite):
        for j in range(Nsite):
            f.write("{:5d} {:5d} {:5d}\n".format(i,j,cnt))
            cnt += 1
    for i in range(Nsite**2):
        f.write("{:5d} {:5d}\n".format(i,1))
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

    file_orbitalidx = "orbitalidx.def"
    make_orbitalidx(dir_def,file_orbitalidx,Nsite)
