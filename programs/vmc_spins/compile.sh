#!/bin/sh

gfortran -o lattice_1d.x lattice_1d.f -llapack -lblas
gfortran -o lattice_2d.x lattice_2d.f -llapack -lblas
gfortran -o preptj_1d.x preptj_1d.f -llapack -lblas
gfortran -o preptj_2d.x preptj_2d.f -llapack -lblas
gcc -c random.c -w   
gfortran -o monte_carlo_tj.x monte_carlo_tj.f random.o -llapack -lblas -fallow-argument-mismatch -w
gfortran -o read_ene.x read_ene.f
gfortran -o read_par.x read_par.f


