#!/usr/bin/bash

N=21
LIST=('1' '2' '4' '8' '16' '10' '20' '5' '11' '13' '17' '19')

for A in "${LIST[@]}"
do
  JULIA_NUM_THREADS=2 julia Shor.jl $N $A
  echo
done
