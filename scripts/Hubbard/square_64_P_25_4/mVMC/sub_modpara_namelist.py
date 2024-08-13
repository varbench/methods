import numpy as np
import json
import os
import time

def make_namelist(output,pre="",post=""):
    ### namelist.def ###
    fileName = []
    fileName.append(pre+'namelist'+post+'.def')
    fileName.append(pre+'modpara'+post+'.def')
    fileName.append(pre+'locspn.def')
    fileName.append(pre+'trans.def')
    fileName.append(pre+'coulombintra.def')
    #fileName.append(pre+'coulombinter.def')
    fileName.append(pre+'greenone.def')
    fileName.append(pre+'greentwo.def')
    #fileName.append(pre+'gutzwiller.def')
    #fileName.append(pre+'jastrowidx.def')
    fileName.append(pre+'orbitalidx.def')
    #fileName.append(pre+'dh4.def')
    fileName.append(pre+'qptransidx.def')
    #fileName.append(pre+'init'+post+'.def')
    fileName.append(pre+'rbm_hlayeridx.def')
    fileName.append(pre+'rbm_playeridx.def')
    fileName.append(pre+'rbm_phidx.def')
    #
    objName = []
    objName.append('         ModPara  ')
    objName.append('         LocSpin  ')
    objName.append('           Trans  ')
    objName.append('    CoulombIntra  ')
    #objName.append('    CoulombInter  ')
    objName.append('       OneBodyG  ')
    objName.append('       TwoBodyG  ')
    #objName.append('      Gutzwiller  ')
    #objName.append('         Jastrow  ')
    objName.append('         Orbital  ')
    #objName.append('             DH4  ')
    objName.append('        TransSym  ')
    #objName.append('       InOrbital  ')
    objName.append('GeneralRBM_HiddenLayer  ')
    objName.append('GeneralRBM_PhysLayer    ')
    objName.append('GeneralRBM_PhysHidden   ')
    #
    f = open(output+"/"+fileName[0],"w")
    for i in range(len(objName)):
        f.write(objName[i]+fileName[i+1]+"\n")
    f.close()

def make_modpara(output,pre="",post="",
                 NVMCCalMode=0,NLanczosMode=0,
                 NDataIdxStart=0,NDataQtySmp=10,
                 Nsite=16,Ne=16,APsgn=0,
                 NSROptItrStep=1000,NSROptItrSmp=100,
                 NVMCWarmUp=100,
                 NVMCInterval=1,
                 NVMCSample=100,
                 NSplitSize=1,rnd=12345,
                 alpha=4):
    ### modpara.def ###
    fileName = []
    fileName.append(pre+'namelist'+post+'.def')
    fileName.append(pre+'modpara'+post+'.def')
    separator = '--------------------\n'
    separator2 = '--------------------'
    f = open(output+"/"+fileName[1],'w')
    f.write(
        separator+
        "Model_Parameters   0\n"+
        separator+
        "VMC_Cal_Parameters\n"+
        separator+
        "CDataFileHead  "+pre+"zvo"+post+"\n"+
        "CParaFileHead  "+pre+"zqp"+post+"\n"+
        separator+
        "NVMCCalMode    {}\n".format(NVMCCalMode)+
        "NLanczosMode   {}\n".format(NLanczosMode)+
        separator+
        "NDataIdxStart  {}\n".format(NDataIdxStart)+
        "NDataQtySmp    {}\n".format(NDataQtySmp)+
        separator+
        "Nsite          {}\n".format(Nsite)+
        "Ncond          {}\n".format(Ne)+
        "2Sz            0\n"+
        "NSPGaussLeg    8\n"+
        "NSPStot        0\n"+
        "NMPTrans       {}\n".format(APsgn*Nsite)+
        "NSROptItrStep  {}\n".format(NSROptItrStep)+
        "NSROptItrSmp   {}\n".format(NSROptItrSmp)+
        "DSROptRedCut   1e-10\n"+
        "DSROptStaDel   1e-3\n"+
        "DSROptStepDt   3e-3\n"+
        "NVMCWarmUp     {}\n".format(NVMCWarmUp)+
        "NVMCInterval   {}\n".format(NVMCInterval)+
        "NVMCSample     {}\n".format(NVMCSample)+
        "NExUpdatePath  0\n"+
#        "NExUpdatePath  1\n"+
        "RndSeed        {}\n".format(rnd)+
        "NSplitSize     {}\n".format(NSplitSize)+
        "NStore         1\n"+
        "NneuronGeneral {}\n".format(alpha*Nsite)
        )
    f.close()

def make_modpara_namelist(dir_def,Nsite,Ne,APsgn,rnd,alpha):
    make_namelist(dir_def,pre="",post="") # optimization
    make_namelist(dir_def,pre="",post="_aft") # calc green func
    make_modpara(dir_def,pre="",post="",NVMCCalMode=0,NLanczosMode=0,NDataIdxStart=1,NDataQtySmp=1,
                 Nsite=Nsite,Ne=Ne,APsgn=APsgn,rnd=rnd,alpha=alpha)
    make_modpara(dir_def,pre="",post="_aft",NVMCCalMode=1,NLanczosMode=1,NDataIdxStart=1,NDataQtySmp=10,
                 Nsite=Nsite,Ne=Ne,APsgn=APsgn,rnd=rnd+10,alpha=alpha,
                 NVMCWarmUp=500,NVMCInterval=5,NVMCSample=5000)
    return 0

if __name__ == "__main__":
    file_json = "dat_output_makedef.json"
    with open(file_json,"r") as json_file:
        loaded_dict = json.load(json_file)
    #print(loaded_dict)

    dir_def = loaded_dict["dir_def"]
    Lx = loaded_dict["Lx"]
    Ly = loaded_dict["Ly"]
    Nsite = Lx*Ly
    Ne = loaded_dict["Ne"]
    APsgn = 1
    if loaded_dict["APFlag"] == 1:
        APsgn = -1;
    alpha = loaded_dict["alpha"]

    current_time = int(time.time())
    process_id = os.getpid()
    random_number = np.random.randint(0, 10000)
    seed = current_time + process_id + random_number
    np.random.seed(seed)
    #rnd = np.random.randint(9999998,size=3)+1 # rnd up to 7 digits and rnd>=1
    rnd = np.random.randint(9999998)+1 # rnd up to 7 digits and rnd>=1
    print("#",current_time,process_id,random_number,seed)
    print("# modpara_rnd",rnd)

    if not os.path.exists(dir_def):
        os.makedirs(dir_def)

    make_modpara_namelist(dir_def,Nsite,Ne,APsgn,rnd,alpha)
