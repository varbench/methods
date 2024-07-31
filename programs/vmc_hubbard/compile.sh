#!/bin/sh

gfortran -o Prep_uniform_wf.x Prep_uniform_wf.f -llapack -lblas
gfortran -o Prep_stripes_wf.x Prep_stripes_wf.f -llapack -lblas
gcc -c Random.c
mpifort -o MC_uniform_wf.x MC_uniform_wf.f Random.o -llapack -lblas -fallow-argument-mismatch
mpifort -o MC_stripe_wf.x MC_stripe_wf.f Random.o -llapack -lblas -fallow-argument-mismatch
gfortran -o Parameters.x Parameters.f
gfortran -o Energy.x Energy.f90
gfortran -o Variance.x Variance.f Random.o
