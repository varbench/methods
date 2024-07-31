#!/bin/sh

cd "$(dirname "$0")"

#Compile codes
cd ../../../programs/vmc_gutzwiller || exit
chmod +x compile.sh
./compile.sh
cd -

#Run VMC
cd vmc_gutzwiller_inputs || exit
../../../../programs/vmc_gutzwiller/lattice_1d.x <wf.in >wf.out
../../../../programs/vmc_gutzwiller/preptj_1d.x <datprep.in >prep.out
../../../../programs/vmc_gutzwiller/monte_carlo_tj.x <VMC.in >VMC.out

#Read data
echo "0 10000 5 0" | ../../../../programs/vmc_gutzwiller/read_ene.x >reading  #energy (<H>/N) + error
echo "<H>/N: $(tail -n 1 fort.20)"

echo "0 10000 5 1" | ../../../../programs/vmc_gutzwiller/read_ene.x >reading #energy**2 (<H**2>/N**2) + error
echo "<H**2>/N**2: $(tail -n 1 fort.20)"

#Note: statistical errors can be reduced by performing more MC steps
