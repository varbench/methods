#!/bin/sh

mkdir -p build
cd build || exit
cmake ..
make
