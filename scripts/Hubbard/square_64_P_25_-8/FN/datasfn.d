&montecarlo                
itestr=1
iseed=19934663
iopt=1                     
iread=0                    
nw=500
ngen=50000                
nscra=150
epst=1.d-10               
/
&fixednode
etry=-0.915d0
tbra=0.08d0
gamma=0.d0
/
&vmc
nbra=64
sflip=0.d0      
/                 
&optimization
nbra=64
sflip=0.d0
nweight=5000
ibinit=50
epsdgel=1.d-06
tau=0.02d0
nraw=0
/
