#!/usr/bin/bash

N=15
LIST=('1' '2' '4' '8' '7' '11' '13')

for A in "${LIST[@]}"
do
  JULIA_NUM_THREADS=2 julia Shor.jl $N $A
  echo
done
