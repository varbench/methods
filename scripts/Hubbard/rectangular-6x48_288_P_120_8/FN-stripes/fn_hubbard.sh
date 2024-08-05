#!/bin/sh

#Compile codes
cd ../../../programs/vmc_hubbard || exit
chmod +x compile.sh
./compile.sh
cd -

#Run VMC
../../../../programs/vmc_hubbard/Prep_stripes_wf.x < datprep.d > prep.out
../../../../programs/vmc_hubbard/Main_stripes_wf.x < datasfn.d > FN.out

#Read data

echo "30 500 2 0" | ../../../../programs/vmc_hubbard/Energy.x  #energy (<H>/N) + error
echo "<H>/N: $(tail -n 1 fort.30)"

#Note: statistical errors can be reduced by performing more MC steps
