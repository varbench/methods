import numpy as np
import json
import os
import time
#
import sub_coulomb
import sub_greenone
import sub_greentwo
import sub_locspn
import sub_modpara_namelist
import sub_orbitalidx
import sub_qptransidx
import sub_rbm_hlayeridx
import sub_rbm_phidx
import sub_rbm_playeridx
import sub_text2dict
import sub_trans

if __name__ == "__main__":
    ## json file
    sub_text2dict.make_text2dict()
    file_json = "dat_output_makedef.json"
    with open(file_json,"r") as json_file:
        loaded_dict = json.load(json_file)
    print(loaded_dict)

    ## load parameters
    dir_def = loaded_dict["dir_def"]
    Lx = loaded_dict["Lx"]
    Ly = loaded_dict["Ly"]
    Nsite = Lx*Ly
    Ne = loaded_dict["Ne"]
    APsgn = 1
    if loaded_dict["APFlag"] == 1:
        APsgn = -1;
    alpha = loaded_dict["alpha"]
    flag_ap = loaded_dict["APFlag"]
    U = loaded_dict["U"]

    ## check directory
    if not os.path.exists(dir_def):
        os.makedirs(dir_def)

    ## output definition files
    file_coulombintra = "coulombintra.def"
    sub_coulomb.make_coulombintra(dir_def,file_coulombintra,Nsite,U)

    file_trans = "trans.def"
    sub_trans.make_trans(dir_def,file_trans,Lx,Ly,Nsite)

    file_greenone = "greenone.def"
    sub_greenone.make_greenone(dir_def,file_greenone,Nsite)

    file_greentwo = "greentwo.def"
    sub_greentwo.make_greentwo(dir_def,file_greentwo,Nsite)

    file_locspn = "locspn.def"
    sub_locspn.make_locspn(dir_def,file_locspn,Nsite)

    file_qptransidx = "qptransidx.def"
    sub_qptransidx.make_qptransidx(dir_def,file_qptransidx,Nsite,Lx,Ly)

    file_rbm = "rbm_hlayeridx.def"
    sub_rbm_hlayeridx.make_rbm(dir_def,file_rbm,Nsite,alpha)

    file_rbm = "rbm_phidx.def"
    sub_rbm_phidx.make_rbm(dir_def,file_rbm,Nsite,Lx,Ly,alpha)

    file_rbm = "rbm_playeridx.def"
    sub_rbm_playeridx.make_rbm(dir_def,file_rbm,Nsite,alpha)

    file_orbitalidx = "orbitalidx.def"
    sub_orbitalidx.make_orbitalidx(dir_def,file_orbitalidx,Nsite)

    current_time = int(time.time())
    process_id = os.getpid()
    random_number = np.random.randint(0, 10000)
    seed = current_time + process_id + random_number
    np.random.seed(seed)
    rnd = np.random.randint(9999998)+1 # rnd up to 7 digits and rnd>=1
    print("#",current_time,process_id,random_number,seed)
    print("# modpara_rnd",rnd)
    sub_modpara_namelist.make_modpara_namelist(dir_def,Nsite,Ne,APsgn,rnd,alpha)
