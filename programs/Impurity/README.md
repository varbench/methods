# Fork Tensor Product States (forkTPS) Library

Fork Tensor Product States is a tensor network states variational wave function tailored for quantum impurity models. This library is used to generate the Ground State Energy and Variance of quantum impurity models in this work.

For detailed information on forkTPS, please refer to [Daniel's paper](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.7.031013).

## Dependencies

This library depends on:
1. [TRIQS](https://triqs.github.io/triqs/latest/install.html) for DMFT-related calculations
2. [ITensor (C++ version)](https://itensor.org/) as the backbone for tensor operations

Both TRIQS and ITensor are open-source libraries. Please follow their respective websites for prerequisite dependencies and installation instructions. Below, we provide the specific versions and installation scripts used in generating data for this work.

## Installation

### 1. TRIQS

Our calculations were performed using TRIQS version 3.1.x. Use the following bash script for installation:


```bash
#!/bin/bash
rm -rf triqs.build
rm -rf triqs.src
# Set this variable to your desired install directory
INSTALL_PREFIX=$(pwd)/install

# Set the number of cores for the compilation
NCORES=16

# Clone the git repository of triqs
git clone -b 3.1.x https://github.com/TRIQS/triqs triqs.src
# Use cmake to configure the triqs build process
mkdir -p triqs.build && cd triqs.build
cmake ../triqs.src -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX

# Build, test and install triqs
make -j$NCORES && make test && make install
cd ../

# Load the triqs installation into your environment
source $INSTALL_PREFIX/share/triqs/triqsvars.sh

echo 
echo "If you want to automatically load triqs into your environment,"
echo "please add the following line to your ~/.bash_profile (or ~/.zprofile):"
echo "source $INSTALL_PREFIX/share/triqs/triqsvars.sh"
```

## ITensor

The C++ version of ITensor was used in this work. The source files are located in the **ITensor_cpp20** folder. To compile the library:

1. Modify `options.mk` according to your specific platform
2. Run the make command:
   ```bash
   make
    ```

## forkTPS

Currently, forkTPS is not publicly available yet. We provide a copy of the library used in this work in the **forkTPS.src** folder. Use the following bash script for installation:


```bash
#!/bin/bash

# Set this variable to your desired install directory
rm -rf forkTPS.build
INSTALL_PREFIX=$(pwd)/install
mkdir -p $INSTALL_PREFIX
LDFLAGS="-L/usr/lib/x86_64-linux-gnu" #set the path to your lapack library
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
```







