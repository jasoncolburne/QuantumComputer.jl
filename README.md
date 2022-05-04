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

<a id='QuantumComputer.jl-Documentation'></a>

<a id='QuantumComputer.jl-Documentation-1'></a>

# QuantumComputer.jl Documentation



<a id='QuantumComputer' href='#QuantumComputer'>#</a>
**`QuantumComputer`** &mdash; *Module*.



```julia
using QuantumComputer
```

a quantum computer simulator. built out of boredom and intrigue. if you are new to julia, be warned - array indices start at 1. throughout this module, qubit index 1 is considered the most significant qubit. in typical circuit diagrams, it represents the top line.

i am fully aware that i 'over-typed' this module. i am new to julia and was curious about performance differences but i assume the compiler will sort everything out in many of the cases i am being explicit.

what this module can do:

  * create registers, superpositions, gates and circuits
  * apply circuits to superpositions
  * add gates to circuits
  * add measurements to circuits
  * add circuits to other circuits as subcircuits
  * add `HybridComponents` to a circuit, which are components where there is a classical interaction with the quantum portion of the simulator
  * measure arbitrary qubits from a superposition
  * generate a variety of single- and multi-qubit static and dynamic gates
  * generate some circuits
  * construct a gate from a circuit
  * retain phase of state of unmeasured qubits (unsure if this is a proper collapse, but was needed for Beauregard's one qubit trick)
  * cache gates on disk

what this module cannot do:

  * decompose a gate into a circuit of smaller gates (planned, long term)
  * print a graphical representation of the circuit (planned, relatively simple but also tedious so not planned for near future)
  * use multiple registers in one circuit (right now you get one quantum and one classical)
  * output bloch sphere animations of algorithms running (planned)

**Example**

```julia
using QuantumComputer

qubit_count = 6
initial_value = 12
constant = 13

register = QuantumComputer.Register(qubit_count, initial_value)
classical_register = QuantumComputer.ClassicalRegister(qubit_count)
superposition = QuantumComputer.Superposition(register.qubits)
adder = QuantumComputer.Circuits.constant_adder(constant, qubit_count)
measurement = QuantumComputer.Measurement(Array(1:qubit_count), Array(1:qubit_count))

circuit = QuantumComputer.Circuit()
QuantumComputer.add_subcircuit_to_circuit!(circuit, adder)
QuantumComputer.add_measurement_to_circuit!(circuit, measurement)

QuantumComputer.apply_circuit_to_superposition!(superposition, circuit, classical_register)

println(classical_register.value)
# 25
```


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L1-L51' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.gate_root_swap' href='#Main.QuantumComputer.gate_root_swap'>#</a>
**`Main.QuantumComputer.gate_root_swap`** &mdash; *Constant*.



[wikipedia](https://en.wikipedia.org/wiki/Quantum_logic_gate#Square_root_of_swap_gate)


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L268-L270' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Circuit' href='#Main.QuantumComputer.Circuit'>#</a>
**`Main.QuantumComputer.Circuit`** &mdash; *Type*.



```julia
Circuit()
```

a quantum circuit. technically this is just an ordered collection of components.


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L671-L675' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.ClassicalRegister' href='#Main.QuantumComputer.ClassicalRegister'>#</a>
**`Main.QuantumComputer.ClassicalRegister`** &mdash; *Type*.



```julia
ClassicalRegister(width, value)
```

a classical register is used to store the results of a measurement

**Arguments**

  * `width`: the width, in bits, of the register
  * `value`: the intitial value of the register (optional, default = 0)


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L61-L69' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Gate' href='#Main.QuantumComputer.Gate'>#</a>
**`Main.QuantumComputer.Gate`** &mdash; *Type*.



```julia
Gate(matrix)
```

a gate operating on a state of size 2^n is really just a 2^n x 2^n unitary matrix. this can be added to a `Circuit`.

**Arguments**

  * `matrix`: a complex unitary matrix of size 2^n, n a positive integer


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L210-L217' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.HybridComponent' href='#Main.QuantumComputer.HybridComponent'>#</a>
**`Main.QuantumComputer.HybridComponent`** &mdash; *Type*.



```julia
HybridComponent(executor, arguments)
```

a component that can be added to a `Circuit` that is capable of controlling quantum gates and subcircuits with classical logic.

**Arguments**

  * `executor`: a function that accepts circuit application arguments (see below for signature) and runs code that applies quantum `Gates` and `Circuits` intelligently
  * `arguments`: fixed arguments known at the time of `Circuit` composition


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L651-L659' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Measurement' href='#Main.QuantumComputer.Measurement'>#</a>
**`Main.QuantumComputer.Measurement`** &mdash; *Type*.



```julia
Measurement(bits_to_output, qubits_to_measure, sample_size)
```

creates a measurement component for incorporation in a `Circuit`. `bits_to_output` and `qubits_to_measure` must be equal in length and contain no duplicates. it is fine for there to be an intersection of the to arrays, as long as each array has unique values.

**Arguments**

  * `bits_to_output`: the bits to write in the classical output register. order must correspond to `qubits_to_measure` order.
  * `qubits_to_measure`: the qubits to measure, in the order of measurement.
  * `sample_size`: the number of samples to perform when measuring


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L490-L499' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Register' href='#Main.QuantumComputer.Register'>#</a>
**`Main.QuantumComputer.Register`** &mdash; *Type*.



```julia
Register(qubit_count, value)
```

a quantum register represents a vector of qubits

**Arguments**

  * `qubit_count`: the width of the register
  * `value`: the initial value of the qubits (optional, default = 0)


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L81-L89' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Superposition' href='#Main.QuantumComputer.Superposition'>#</a>
**`Main.QuantumComputer.Superposition`** &mdash; *Type*.



```julia
Superposition(qubits)
Superposition(state)
```

a superposition is a representation of the probabilities of all possible states of a vector of qubits, typically a register (but there is no hard requirement for a register as input in this simulator)


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L109-L114' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.add_gate_to_circuit!-Tuple{Main.QuantumComputer.Circuit,Main.QuantumComputer.Gate}' href='#Main.QuantumComputer.add_gate_to_circuit!-Tuple{Main.QuantumComputer.Circuit,Main.QuantumComputer.Gate}'>#</a>
**`Main.QuantumComputer.add_gate_to_circuit!`** &mdash; *Method*.



```julia
add_gate_to_circuit!(circuit, gate)
```

adds a `Gate` to a `Circuit`.

**Arguments**

  * `circuit`: the `Circuit`
  * `gate`: the `Gate`


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L684-L692' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.add_hybrid_component_to_circuit!-Tuple{Main.QuantumComputer.Circuit,Main.QuantumComputer.HybridComponent}' href='#Main.QuantumComputer.add_hybrid_component_to_circuit!-Tuple{Main.QuantumComputer.Circuit,Main.QuantumComputer.HybridComponent}'>#</a>
**`Main.QuantumComputer.add_hybrid_component_to_circuit!`** &mdash; *Method*.



```julia
add_hybrid_component_to_circuit!(circuit, component)
```

adds a `HybridComponent` to a `Circuit`.

**Arguments**

  * `circuit`: the `Circuit`
  * `component`: the `HybridComponent`


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L710-L718' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.add_measurement_to_circuit!-Tuple{Main.QuantumComputer.Circuit,Main.QuantumComputer.Measurement}' href='#Main.QuantumComputer.add_measurement_to_circuit!-Tuple{Main.QuantumComputer.Circuit,Main.QuantumComputer.Measurement}'>#</a>
**`Main.QuantumComputer.add_measurement_to_circuit!`** &mdash; *Method*.



```julia
add_measurement_to_circuit!(circuit, measurement)
```

adds a `Measurement` to a `Circuit`.

**Arguments**

  * `circuit`: the `Circuit`
  * `measurement`: the `Measurement`


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L697-L705' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.add_subcircuit_to_circuit!-Tuple{Main.QuantumComputer.Circuit,Main.QuantumComputer.Circuit}' href='#Main.QuantumComputer.add_subcircuit_to_circuit!-Tuple{Main.QuantumComputer.Circuit,Main.QuantumComputer.Circuit}'>#</a>
**`Main.QuantumComputer.add_subcircuit_to_circuit!`** &mdash; *Method*.



```julia
add_subcircuit_to_circuit!(circuit, subcircuit)
```

adds a `Circuit` to a another `Circuit` as a component.

**Arguments**

  * `circuit`: the parent `Circuit`
  * `subcircuit`: the child `Circuit`


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L723-L731' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.apply_circuit_to_superposition!' href='#Main.QuantumComputer.apply_circuit_to_superposition!'>#</a>
**`Main.QuantumComputer.apply_circuit_to_superposition!`** &mdash; *Function*.



```julia
apply_circuit_to_superposition!(superposition, circuit, classical_register)
```

applies the gates and measurements in a circuit to the superposition provided. a classical register is required for measurement outputs.

**Arguments**

  * `superposition`: the `Superposition` to operate on
  * `circuit`: the `Circuit` to apply
  * `classical_register`: a `ClassicalRegister` capable of holding all circuit measurement output bits


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L820-L829' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.gate_cnx-Tuple{Array{Int64,N} where N,Int64,Int64}' href='#Main.QuantumComputer.gate_cnx-Tuple{Array{Int64,N} where N,Int64,Int64}'>#</a>
**`Main.QuantumComputer.gate_cnx`** &mdash; *Method*.



```julia
gate_cnx(control_qubits, qubit_index, qubit_count)
```

n-controlled pauli x gate.

**Arguments**

  * `control_qubits`: the indices of the qubits controlling the pauli x operation
  * `qubit_index`: the index of the qubit upon which to perform the pauli x operation
  * `qubit_count`: the total number of qubits the resultant gate acts on


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L418-L427' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.gate_control-Tuple{Main.QuantumComputer.Gate,Array{Int64,N} where N,Int64,Int64}' href='#Main.QuantumComputer.gate_control-Tuple{Main.QuantumComputer.Gate,Array{Int64,N} where N,Int64,Int64}'>#</a>
**`Main.QuantumComputer.gate_control`** &mdash; *Method*.



```julia
gate_control(gate, control_qubits, qubit_index, qubit_count)
```

creates an n-controlled single qubit gate. acts on an arbitrary qubit in a superposition and is controlled by n arbitrary qubits.

**Arguments**

  * `gate`: the underlying single qubit gate
  * `control_qubits`: an array of integers selecting control qubits
  * `qubit_index`: an integer selecting the qubit to which the single-qubit gate is applied
  * `qubit_count`: the number of qubits in the superposition this gate will act on


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L354-L364' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.gate_cx-Tuple{Int64,Int64,Int64}' href='#Main.QuantumComputer.gate_cx-Tuple{Int64,Int64,Int64}'>#</a>
**`Main.QuantumComputer.gate_cx`** &mdash; *Method*.



```julia
gate_cx(control_qubit, qubit_index, qubit_count)
```

controlled pauli x gate.

**Arguments**

  * `control_qubit`: the index of the qubit controlling the pauli x operation
  * `qubit_index`: the index of the qubit upon which to perform the pauli x operation
  * `qubit_count`: the total number of qubits the resultant gate acts on


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L404-L413' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.gate_extension-Tuple{Main.QuantumComputer.Gate,Int64,Int64}' href='#Main.QuantumComputer.gate_extension-Tuple{Main.QuantumComputer.Gate,Int64,Int64}'>#</a>
**`Main.QuantumComputer.gate_extension`** &mdash; *Method*.



```julia
gate_extension(gate, qubit_index, qubit_count)
```

a gate that acts on a subset of qubits by identifying a range of qubits in a superposition and applying a smaller gate to that range, while applying identity transformations to the qubits outside the range. effectively, extends a gate smaller than required by the  superposition size to operate on a subset of qubits in that superposition.

**Arguments**

  * `gate`: the underlying gate
  * `qubit_index`: an integer selecting the qubit where the underlying gate begins
  * `qubit_count`: the number of qubits in the superposition this gate will act on


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L329-L342' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.gate_fourier_transform-Tuple{Int64}' href='#Main.QuantumComputer.gate_fourier_transform-Tuple{Int64}'>#</a>
**`Main.QuantumComputer.gate_fourier_transform`** &mdash; *Method*.



```julia
gate_fourier_transform(qubit_count)
```

the quantum fourier transform. [wikipedia](https://en.wikipedia.org/wiki/Quantum_Fourier_transform)

**Arguments**

  * `qubit_count`: the number of qubits the gate operates on


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L453-L460' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.gate_invert-Tuple{Main.QuantumComputer.Gate}' href='#Main.QuantumComputer.gate_invert-Tuple{Main.QuantumComputer.Gate}'>#</a>
**`Main.QuantumComputer.gate_invert`** &mdash; *Method*.



```julia
gate_invert(gate)
```

inverts a unitary matrix quickly, taking advantage of the fact that the conjugate transpose is the inverse of a unitary matrix

**Arguments**

  * `gate`: the gate to invert


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L478-L485' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.gate_multi_control-Tuple{Main.QuantumComputer.Gate,Int64,Int64}' href='#Main.QuantumComputer.gate_multi_control-Tuple{Main.QuantumComputer.Gate,Int64,Int64}'>#</a>
**`Main.QuantumComputer.gate_multi_control`** &mdash; *Method*.



```julia
gate_multi_control(gate, control_count, qubit_count)
```

this special function creates a controlled-n gate efficiently, by utilizing the top n qubits

**Arguments**

  * `gate`


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L432-L439' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.gate_swap-Tuple{Int64,Int64,Int64}' href='#Main.QuantumComputer.gate_swap-Tuple{Int64,Int64,Int64}'>#</a>
**`Main.QuantumComputer.gate_swap`** &mdash; *Method*.



```julia
gate_swap(qubit_a_index, qubit_b_index, qubit_count)
```

a gate that swaps two arbitrary qubits in a superposition [explanation](https://quantumcomputing.stackexchange.com/questions/9181/swap-gate-on-2-qubits-in-3-entangled-qubit-system)

**Arguments**

  * `qubit_a_index`: the index of the first qubit to swap
  * `qubit_b_index`: the index of the second qubit to swap
  * `qubit_count`: the total number of qubits in the superposition

this gate will operate on


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L274-L284' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.ket_bra-Tuple{Int64,Int64}' href='#Main.QuantumComputer.ket_bra-Tuple{Int64,Int64}'>#</a>
**`Main.QuantumComputer.ket_bra`** &mdash; *Method*.



```julia
ket_bra(i, j)
```

check the stackexchange post on the previous function to understand what is happening here.

**Arguments**

  * `i`: first value
  * `j`: second value


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L311-L319' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.measure_superposition-Tuple{Main.QuantumComputer.Superposition,Main.QuantumComputer.ClassicalRegister,Main.QuantumComputer.Measurement}' href='#Main.QuantumComputer.measure_superposition-Tuple{Main.QuantumComputer.Superposition,Main.QuantumComputer.ClassicalRegister,Main.QuantumComputer.Measurement}'>#</a>
**`Main.QuantumComputer.measure_superposition`** &mdash; *Method*.



```julia
measure_superposition(superposition, classical_register, measurement)
```

measure some qubits and store the result in a classical register

**Arguments**

  * `superposition`: the superposition being measured
  * `classical_register`: the output register
  * `measurement`: details of the measurement


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L543-L552' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.qubit_tensor_product-Tuple{Array{Complex{Float64},2}}' href='#Main.QuantumComputer.qubit_tensor_product-Tuple{Array{Complex{Float64},2}}'>#</a>
**`Main.QuantumComputer.qubit_tensor_product`** &mdash; *Method*.



```julia
qubit_tensor_product(qubits)
```

this custom function computes a tensor product on the n rows of a matrix of 2-column qubits. this function is threaded, and the core of the algorithm follows in the next function. this works by iteratively applying the tensor product qubit by qubit.

**Arguments**

  * `qubits`: a matrix of qubits, in the same arrangement as a `Register`


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L127-L134' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.qubit_tensor_product_thread-Tuple{Array{Complex{Float64},1},Array{Complex{Float64},1},Array{Complex{Float64},1},Int64,Int64,Int64}' href='#Main.QuantumComputer.qubit_tensor_product_thread-Tuple{Array{Complex{Float64},1},Array{Complex{Float64},1},Array{Complex{Float64},1},Int64,Int64,Int64}'>#</a>
**`Main.QuantumComputer.qubit_tensor_product_thread`** &mdash; *Method*.



```julia
qubit_tensor_product_thread(next_product, basis, product, product_size, thread_number, thread_count)
```

the threaded portion from our qubit tensor product algorithm. this will compute a slice of the current tensor product expansion from the parent algorithm. this fills in a portion of `next_product` by computing a partial tensor product of `basis` and `product`.

**Arguments**

  * `next_product`: the product (state) that is currently being computed
  * `basis`: the qubit being applied to the previous product
  * `product`: the previous product
  * `product_size`: the number of states represented by `product`
  * `thread_number`: the algorithm-assigned thread number of this thread
  * `thread_count`: the total number of threads


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L169-L182' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.qubits_operated_on_by_unitary_matrix-Tuple{Array{T,2} where T}' href='#Main.QuantumComputer.qubits_operated_on_by_unitary_matrix-Tuple{Array{T,2} where T}'>#</a>
**`Main.QuantumComputer.qubits_operated_on_by_unitary_matrix`** &mdash; *Method*.



```julia
qubits_operated_on_by_unitary_matrix(matrix)
```

a helper function to compute the number of qubits a gate operates on

**Arguments**

  * `matrix`: a unitary matrix


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L227-L234' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Circuits' href='#Main.QuantumComputer.Circuits'>#</a>
**`Main.QuantumComputer.Circuits`** &mdash; *Module*.



```julia
QuantumComputer.Circuits
```

a collection of simple quantum circuits.

**Example**

```julia
using QuantumComputer

qubit_count = 6
initial_value = 12
constant = 13

register = QuantumComputer.Register(qubit_count, initial_value)
classical_register = QuantumComputer.ClassicalRegister(qubit_count)
superposition = QuantumComputer.Superposition(register.qubits)
adder = QuantumComputer.Circuits.constant_adder(constant, qubit_count)
measurement = QuantumComputer.Measurement(Array(1:qubit_count), Array(1:qubit_count))

circuit = QuantumComputer.Circuit()
QuantumComputer.add_subcircuit_to_circuit!(circuit, adder)
QuantumComputer.add_measurement_to_circuit!(circuit, measurement)

QuantumComputer.apply_circuit_to_superposition!(superposition, circuit, classical_register)

println(classical_register.value)
# 25
```


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L850-L878' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Circuits.constant_adder-Tuple{Int64,Int64}' href='#Main.QuantumComputer.Circuits.constant_adder-Tuple{Int64,Int64}'>#</a>
**`Main.QuantumComputer.Circuits.constant_adder`** &mdash; *Method*.



```julia
constant_adder(n, qubit_count)
```

a quantum circuit that adds `n` to the superposition's value (`mod 2^qubit_count`)

**Arguments:**

  * `n`: the constant to add
  * `qubit_count`: the number of qubits in the superposition


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L884-L892' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Circuits.constant_adder_core-Tuple{Int64,Int64}' href='#Main.QuantumComputer.Circuits.constant_adder_core-Tuple{Int64,Int64}'>#</a>
**`Main.QuantumComputer.Circuits.constant_adder_core`** &mdash; *Method*.



```julia
constant_adder_core(n, qubit_count)
```

a quantum circuit that adds `n` to the superposition's value (`mod 2^qubit_count`) without applying pre and post fourier transforms

**Arguments:**

  * `n`: the constant to add
  * `qubit_count`: the number of qubits in the superposition


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L907-L915' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Circuits.period_finding_for_11x_mod_15-Tuple{}' href='#Main.QuantumComputer.Circuits.period_finding_for_11x_mod_15-Tuple{}'>#</a>
**`Main.QuantumComputer.Circuits.period_finding_for_11x_mod_15`** &mdash; *Method*.



```julia
period_finding_for_11x_mod_15()
```

a 5 qubit circuit that returns a value related to the period of f(x) = 11^x mod 15. this effectively allows implementation of Shor's algorithm for this specific case.


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L1220-L1224' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Circuits.shor2n3_controlled_controlled_modular_adder-Tuple{Int64,Int64}' href='#Main.QuantumComputer.Circuits.shor2n3_controlled_controlled_modular_adder-Tuple{Int64,Int64}'>#</a>
**`Main.QuantumComputer.Circuits.shor2n3_controlled_controlled_modular_adder`** &mdash; *Method*.



```julia
shor2n3_controlled_controlled_modular_adder(n, a, qubit_count)
```

a quantum circuit that adds `a` to the superposition's value (`mod n`). see the modular adder circuit in [this paper](https://arxiv.org/pdf/quant-ph/0205095.pdf).

**Arguments:**

  * `n`: the modulus
  * `a`: the constant to add


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L936-L944' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Circuits.shor2n3_controlled_modular_multiplier-Tuple{Int64,Int64}' href='#Main.QuantumComputer.Circuits.shor2n3_controlled_modular_multiplier-Tuple{Int64,Int64}'>#</a>
**`Main.QuantumComputer.Circuits.shor2n3_controlled_modular_multiplier`** &mdash; *Method*.



```julia
shor2n3_controlled_modular_multiplier(n, a)
```

see the modular multiplier circuit in [this paper](https://arxiv.org/pdf/quant-ph/0205095.pdf).

**Arguments:**

  * `n`: the modulus
  * `a`: the constant to multiply by


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L1022-L1030' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Circuits.shor2n3_controlled_ua-Tuple{Int64,Int64}' href='#Main.QuantumComputer.Circuits.shor2n3_controlled_ua-Tuple{Int64,Int64}'>#</a>
**`Main.QuantumComputer.Circuits.shor2n3_controlled_ua`** &mdash; *Method*.



```julia
shor2n3_controlled_ua(n, a)
```

the controlled-Ua gate from [this paper](https://arxiv.org/pdf/quant-ph/0205095.pdf).

**Arguments**

  * `n`: the modulus
  * `a`: the constant to multiply `x` by


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L1072-L1080' class='documenter-source'>source</a><br>

<a id='Main.QuantumComputer.Circuits.shor2n3_period_finding-Tuple{Int64,Int64}' href='#Main.QuantumComputer.Circuits.shor2n3_period_finding-Tuple{Int64,Int64}'>#</a>
**`Main.QuantumComputer.Circuits.shor2n3_period_finding`** &mdash; *Method*.



```julia
shor2n3_period_finding(n::Int64, a::Int64)
```

Beauregard's circuit for finding the period of a^x mod n.

**Arguments**

  * `n`: the modulus
  * `a`: the base


<a target='_blank' href='https://github.com/jasoncolburne/QuantumComputer.jl/blob/a0110d431a2ecb6d3be8fabfa93a8648bcef64c6/src/QuantumComputer.jl#L1168-L1176' class='documenter-source'>source</a><br>

