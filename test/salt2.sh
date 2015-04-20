#!/bin/bash

dir=../bin

# salt bridges between residues when pos/neg distance is less than 4.5 Angs. 
# for at least 0.1 (10%) of the trajectory frames
$dir/eucb -mol protein -salt2 -cutoff 0.1,4.5 -seq A,C -last 100

