#!/bin/bash
gfortran -O3 -fcheck=all -flto -c BUF.f90
gfortran -O3 -fcheck=all -flto -c constants.f90
gfortran -O3 -fcheck=all -flto -c sto.f90 
gfortran -O3 -fcheck=all -flto -c linal.f90 
gfortran -O3 -fcheck=all -flto -c orb.f90 linal.o 
gfortran -O3 -fcheck=all xvmc.f90 *.o -o /Users/james.thorpe/work/QMC/bin/xvmc -lblas -llapack
