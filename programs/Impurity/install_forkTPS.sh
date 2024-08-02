#!/bin/bash

# Set this variable to your desired install directory
rm -rf forkTPS.build
INSTALL_PREFIX=$(pwd)/install
mkdir -p $INSTALL_PREFIX
LDFLAGS="-L/usr/lib/x86_64-linux-gnu" #set the path to the lapack library
# Set the number of cores for the compilation
NCORES=16

# Set the path to the ITensor library
ITENSOR_LIBRARIES=/MyLibs/forkTPSLibs/ITensor_cpp20/lib/libitensor-g.a
ITENSOR_INCLUDE_DIR=/MyLibs/forkTPSLibs/ITensor_cpp20/

# Clone the git repository of forkTPS finiteTDLR branch, currently forkTPS is not publickly available yet, we provide a copy of the library in the folder forkTPS.src
# git clone -b finiteTDLR  https://github.com/TRIQS/forktps.git forkTPS.src

# Use cmake to configure the triqs build process
mkdir -p forkTPS.build && cd forkTPS.build

cmake ../forkTPS.src -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX -DCMAKE_LD_FLAGS=-L/usr/lib/x86_64-linux-gnu -DITENSOR_LIBRARIES=$ITENSOR_LIBRARIES -DITENSOR_INCLUDE_DIR=$ITENSOR_INCLUDE_DIR

# Build, test and install triqs
make -j$NCORES && make test && make install
cd ../

