#!/bin/sh

gfortran -o Prep_uniform_wf.x Prep_uniform_wf.f -llapack -lblas
gfortran -o Prep_stripes_wf.x Prep_stripes_wf.f -llapack -lblas
gcc -c Random.c
gfortran -o Main_uniform_wf.x Main_uniform_wf.f Random.o -llapack -lblas -fallow-argument-mismatch
gfortran -o Main_stripes_wf.x Main_stripes_wf.f Random.o -llapack -lblas -fallow-argument-mismatch
gfortran -o Parameters.x Parameters.f
gfortran -o Energy.x Energy.f90
gfortran -o Variance.x Variance.f Random.o
