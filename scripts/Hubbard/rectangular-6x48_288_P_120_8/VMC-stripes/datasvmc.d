&montecarlo                 
itestr=2
iseed=19934663
iopt=1                      
iread=0                     
nw=1
ngen=500000                 
nscra=200
epst=1.d-10                
/
&fixednode
etry=0.d0
tbra=0.1d0
gamma=0.d0
/
&vmc
nbra=288
sflip=0.d0      
/                 
&optimization
nbra=288
sflip=0.d0
nweight=5000
ibinit=50
epsdgel=1.d-06
tau=0.02d0
nraw=0
/
