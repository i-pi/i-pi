#!/bin/sh

f2py --fcompiler=gfortran --f90flags=-ffree-line-length-0 -m fortran -c fortran.f90 
