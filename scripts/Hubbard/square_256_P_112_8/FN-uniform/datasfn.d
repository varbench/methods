&montecarlo                
itestr=1
iseed=19934663
iopt=1                     
iread=0                   
nw=500
ngen=15000                
nscra=200
epst=1.d-10                 
/
&fixednode
etry=-0.755d0
tbra=0.08d0
gamma=0.d0
/
&vmc
nbra=256
sflip=0.2d0     
/                 
&optimization
nbra=256
sflip=0.d0
nweight=312
ibinit=10
epsdgel=1.d-06
tau=0.02d0
nraw=0
/
