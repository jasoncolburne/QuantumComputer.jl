#!/usr/bin/bash

LIST=('2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12')

for N in "${LIST[@]}"
do
  JULIA_NUM_THREADS=2 julia Grover.jl $N
  echo
done
