import numpy as np
import json
import os

def make_rbm(dir_def,file_rbm,Nsite,alpha):
    f = open(dir_def+"/"+file_rbm,"w")
    f.write("====\n")
    f.write("NRBM_HiddenLayerIdx "+"{}\n".format(alpha))
    f.write("ComplexType 1\n")
    f.write("k RBM_HiddenLayer_Idx\n")
    f.write("====\n")
    cnt = 0
    for i in range(alpha):
        for j in range(Nsite):
            f.write("{:d} {:d}\n".format(cnt,i))
            cnt += 1
    for i in range(alpha):
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

    file_rbm = "rbm_hlayeridx.def"
    make_rbm(dir_def,file_rbm,Nsite,alpha)
