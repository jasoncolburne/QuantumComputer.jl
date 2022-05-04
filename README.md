# Getting Started

## Running Tests
```
QuantumComputer.jl$ julia
julia> ]
pkg> activate .
pkg> test
```

## Examples

The `shor15` and `grover` example scripts should run well on a modern laptop with 16gb of RAM (if `grover` won't run to completion, remove the last element - 12 - from the list). To get `shor21` to work on a 16gb computer, I had to pre-build all the required gates in a very specific order and cache them to disk, programatically, one at a time. I suppose I could write a script to do this so that others wishing to do the same on limited hardware can experience the glory of Beauregard's circuit on 13 qubits. I was running a WSL VM inside Windows, so if you were to put something like Debian slim on bare metal you could probably build shor21 without pre-computation with 16gb.

```
examples$ ./grover.sh | tee my.log
examples$ cat my.log
examples$ diff -u my.log grover.txt

examples$ ./shor15.sh # etc
```

[Documentation](DOCUMENTATION.md) for the source code.
