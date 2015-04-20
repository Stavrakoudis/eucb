#!/bin/bash

dir=../bin

# all torsion angles of chain A residues, results stored at rmsf_A.dat rmsf_C.dat files
$dir/eucb -mol protein -tors all -seq A -last 100

# phi,psi torsion angles of chain A arginine residues, uncomment to run
#$dir/eucb -mol protein -tors phi,psi -seq A -resname ARG -last 100
