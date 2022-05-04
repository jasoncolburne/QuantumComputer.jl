#!/usr/bin/bash

julia -e 'using Pkg; Pkg.activate("."); include("src/QuantumComputer.jl"); include("docs/make.jl")'
cp docs/build/index.md DOCUMENTATION.md
