#!/bin/sh

f2py -c --f90exec=mpif90 --f77exec=mpif90 -m vort_div vort_div_openmpi.F90 
