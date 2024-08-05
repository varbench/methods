#!/bin/sh

cd "$(dirname "$0")"

#Compile codes
cd ../../../programs/vmc_hubbard || exit
chmod +x compile.sh
./compile.sh
cd -

#Run VMC
cd vmc_hubbard_inputs || exit
../../../../programs/vmc_hubbard/Prep_stripes.x <datprep.d >prep.out
../../../../programs/vmc_hubbard/Main_stripes.x <datasvmc.d >VMC.out

#Read data
echo "0 1000 5 1" | ../../../../programs/vmc_hubbard/Energy.x  #energy (<H>/N) + error
echo "<H>/N: $(tail -n 1 fort.30)"

echo "0 1000 5 1" | ../../../../programs/vmc_hubbard/Energy.x >#energy**2 (<H**2>/N**2) + error
echo "<H**2>/N**2: $(tail -n 1 fort.30)"

#Note: statistical errors can be reduced by performing more MC steps
