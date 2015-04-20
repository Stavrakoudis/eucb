#!/bin/bash

dir=../bin

# RMSF of CA atoms, results stored at rmsf_A.dat rmsf_C.dat files
$dir/eucb -mol protein -rmsf -last 100

# RMSF of CB atoms if chain A, uncomment the line
#$dir/eucb -mol protein -rmsf -seq A -atomname CB -last 100

# RMSD of backbone atoms, results stored at rmsd.dat
$dir/eucb -mol protein -rmsd -atomname N,CA,C -seq A -last 100

# RMSD of backbone 4 atoms of chain A, uncomment to run
#$dir/eucb -mol protein -rmsd -atomname N,CA,C,O -seq A -last 100
