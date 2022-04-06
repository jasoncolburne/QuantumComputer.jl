"""
    using QuantumComputer

a quantum computer simulator. built out of boredom and intrigue.
if you are new to julia, be warned - array indices start at 1.
throughout this module, qubit index 1 is considered the most
significant qubit. in typical circuit diagrams, it represents the
top line.

what this module can do:
- create registers, superpositions, gates and circuits
- apply circuits to superpositions
- add gates to circuits
- add measurements to circuits
- add circuits to other circuits as subcircuits
- add `HybridComponents` to a circuit, which are components where there
is a classical interaction with the quantum portion of the simulator
- measure arbitrary qubits from a superposition (collapse implementation
will currently lose information for entangled and unmeasured qubits)
- generate a variety of single- and multi-qubit static and dynamic gates
- generate some basic circuits

what this module cannot do:
- properly collapse a superposition during measurement (the squaring/
rooting when we sum the unrealized state into the realized state causes
a loss of information) (need to sum in a more sophistocated way?
is this possible?)
- construct a gate from a circuit (planned, pretty simple to implement)
- decompose a gate into a circuit of smaller gates (unplanned)
- print a graphical representation of the circuit (planned, relatively
simple but also tedious so not planned for near future)
- use multiple registers in one circuit (right now you get one quantum
and one classical)

# Example
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
"""
module QuantumComputer

using LinearAlgebra
using StatsBase

# base components (registers and superposition)
###############################################

"""
    ClassicalRegister(width, value)

a classical register is used to store the results of a measurement

# Arguments
- `width`: the width, in bits, of the register
- `value`: the intitial value of the register (optional, default = 0)
"""
mutable struct ClassicalRegister
  width::Int64
  value::Int64

  function ClassicalRegister(width::Int64, value::Int64 = 0)
    value < 2^width || throw(DomainError(value, "value does not fit in register"))

    new(width, value)
  end
end

"""
    Register(qubit_count, value)

a quantum register represents a vector of qubits

# Arguments
- `qubit_count`: the width of the register
- `value`: the initial value of the qubits (optional, default = 0)
"""
mutable struct Register
  qubits::Matrix{Complex{Float64}}

  function Register(qubit_count::Int64, value::Int64 = 0)
    value < 2^qubit_count || throw(DomainError(value, "value does not fit in register"))

    qubits::Matrix{Int64} = hcat(ones(Int64, qubit_count), zeros(Int64, qubit_count))
    
    for i::Int64 in 1:qubit_count
      if (value & 2^(i - 1)) != 0
        k::Int64 = (qubit_count - i + 1)
        qubits[k:k, :] = [0 1]
      end  
    end

    new(qubits)
  end
end

"""
    Superposition(qubits)
    Superposition(state)

a superposition is a representation of the probabilities of all possible
states of a vector of qubits, typically a register (but there is no hard
requirement for a register as input in this simulator)
"""
mutable struct Superposition
  state::Array{Complex{Float64}, 1}

  function Superposition(qubits::Matrix{Complex{Float64}})
    new(qubit_tensor_product(qubits))
  end

  function Superposition(state::Array{Complex{Float64}, 1})
    new(state)
  end
end

"""
    qubit_tensor_product(qubits)

this custom function computes a tensor product on the n rows of
a matrix of 2-column qubits. this function is threaded, and the
core of the algorithm follows in the next function. this works by
iteratively applying the tensor product qubit by qubit.

# Arguments
- `qubits`: a matrix of qubits, in the same arrangement as a `Register`
"""
function qubit_tensor_product(qubits::Matrix{Complex{Float64}})::Array{Complex{Float64}, 1}
  qubit_count::Int64 = size(qubits, 1)
  thread_count::Int64 = Threads.nthreads()

  product::Array{Complex{Float64}, 1} = qubit_count > 1 ? qubits[1, :] : throw(DomainError(qubits, "at least 2 qubits required"))
  for i::Int64 in 2:qubit_count
    basis::Array{Complex{Float64}, 1} = qubits[i, :]
    product_size::Int64 = length(product)
    next_product::Array{Complex{Float64}, 1} = Array{Complex{Float64}}(undef, 2 * product_size)
    
    # coming back to this code, i was confused why we're iterating to 
    # min(thread_count, product_size) but passing the child function thread_count
    # turns out that since we need to do some computations with product size in
    # the child function, we end up handling this there.
    Threads.@threads for thread_number::Int64 in 1:min(thread_count, product_size)
      qubit_tensor_product_thread(next_product, basis, product, product_size, thread_number, thread_count)
    end
    
    product = next_product
  end

  product
end

"""
    qubit_tensor_product_thread(next_product, basis, product, product_size, thread_number, thread_count)

the threaded portion from our qubit tensor product algorithm. this
will compute a slice of the current tensor product expansion from the
parent algorithm.

this fills in a portion of `next_product` by computing a partial
tensor product of `basis` and `product`.

# Arguments
- `next_product`: the product (state) that is currently being computed
- `basis`: the qubit being applied to the previous product
- `product`: the previous product
- `product_size`: the number of states represented by `product`
- `thread_number`: the algorithm-assigned thread number of this thread
- `thread_count`: the total number of threads
"""
function qubit_tensor_product_thread(next_product::Array{Complex{Float64}, 1}, basis::Array{Complex{Float64}, 1}, product::Array{Complex{Float64}, 1}, product_size::Int64, thread_number::Int64, thread_count::Int64)
  width::Int64, padded_count::Int64 = divrem(product_size, thread_count)

  if padded_count < thread_number
    initial_value = (width + 1) * padded_count + width * (thread_number - padded_count - 1) + 1
  else
    width += 1
    initial_value = width * (thread_number - 1) + 1
  end

  for i::Int64 in initial_value:min(initial_value + width - 1, product_size)
    next_product[2i - 1] = basis[1] * product[i]
    next_product[2i] = basis[2] * product[i]
  end
end

# gates
#######

"""
    Gate(matrix)

a gate operating on a state of size 2^n is really just a 2^n x 2^n unitary matrix.
this can be added to a `Circuit`.

# Arguments
- `matrix`: a complex unitary matrix of size 2^n, n a positive integer
"""
struct Gate
  matrix::Matrix{Complex{Float64}}
  superposed_qubits_required::Int64

  function Gate(matrix::Matrix)
    new(matrix, qubits_operated_on_by_unitary_matrix(matrix))
  end
end

"""
    qubits_operated_on_by_unitary_matrix(matrix)

a helper function to compute the number of qubits a gate operates on

# Arguments
- `matrix`: a unitary matrix
"""
function qubits_operated_on_by_unitary_matrix(matrix::Matrix)::Int64
  dimensions = size(matrix)
  dimensions[1] == dimensions[2] || throw(DomainError(matrix, "matrix must be square"))
  log2(dimensions[1])
end

# single qubit gates
####################

# https://qiskit.org/textbook/ch-states/single-qubit-gates.html

# pauli
gate_x = Gate([0 1; 1 0])
gate_y = Gate([0.0 -1.0im; 1.0im 0.0])
gate_z = Gate([1 0; 0 -1])
# hadamard
gate_h = Gate([1 1; 1 -1] / √2)
# phase
gate_p(ϕ::Float64) = Gate([1 0; 0 ℯ^(ϕ*im)])
# identity
gate_i = Gate([1 0; 0 1])
# pi / 2 and pi / 4 phase rotations
gate_s = gate_p(pi/2)
gate_sdg = gate_p(-pi/2)
gate_t = gate_p(pi/4)
gate_tdg = gate_p(-pi/4)
# u, the most general gate
gate_u(θ::Float64, ϕ::Float64, λ::Float64) = Gate([cos(θ/2) -sin(θ/2)*ℯ^(λ*im); sin(θ/2)*ℯ^(ϕ*im) cos(θ/2)*ℯ^((ϕ+λ)*im)])

# multi-qubit gates
###################

"""
https://en.wikipedia.org/wiki/Quantum_logic_gate#Square_root_of_swap_gate
"""
gate_root_swap = Gate([1 0 0 0; 0 (1.0+im)/2 (1.0-im)/2 0; 0 (1.0-im)/2 (1.0+im)/2 0; 0 0 0 1])

"""
    gate_swap(qubit_a_index, qubit_b_index, qubit_count)

https://quantumcomputing.stackexchange.com/questions/9181/swap-gate-on-2-qubits-in-3-entangled-qubit-system
a gate that swaps two arbitrary qubits in a superposition

# Arguments
- `qubit_a_index`: the index of the first qubit to swap
- `qubit_b_index`: the index of the second qubit to swap
- `qubit_count`: the total number of qubits in the superposition
this gate will operate on
"""
function gate_swap(qubit_a_index::Int64, qubit_b_index::Int64, qubit_count::Int64)
  identity = ((1.0 + 0.0im)*I)(2)

  matrix_size = 2^qubit_count
  matrix::Matrix{Complex{Float64}} = zeros(Complex{Float64}, (matrix_size, matrix_size))
  matrices::Array{Matrix{Complex{Float64}}, 1} = Array{Matrix{Complex{Float64}}}(undef, 0)
  
  for i in 1:qubit_count
    push!(matrices, identity)
  end

  for i in 0:1
    for j in 0:1
      matrices[qubit_a_index] = ket_bra(i, j)
      matrices[qubit_b_index] = ket_bra(j, i)
      partial_sum = matrices[1]
      for k in 2:qubit_count
        partial_sum = kron(partial_sum, matrices[k])
      end
      matrix += partial_sum
    end
  end

  Gate(matrix)
end

"""
    ket_bra(i, j)

check the stackexchange post on the previous function to understand
what is happening here.

# Arguments
- `i`: first value
- `j`: second value
"""
function ket_bra(i::Int64, j::Int64)::Matrix{Complex{Float64}}
  (i > 1 || i < 0 || j > 1 || j < 0) && throw(DomainError((i, j), "arguments must be in range [0, 1]"))

  matrix::Matrix{Complex{Float64}} = zeros(Complex{Float64}, (2, 2))
  matrix[[(j << 1) + i + 1]] = [1.0 + 0.0im]
  matrix
end

"""
    gate_extension()

a gate that acts on a subset of qubits by identifying a range of
qubits in a superposition and applying a smaller gate to that range,
while applying identity transformations to the qubits outside the
range. effectively, extends a gate smaller than required by the 
superposition size to operate on a subset of qubits in that superposition.

# Arguments
- `gate`: the underlying gate
- `qubit_index`: an integer selecting the qubit where the underlying gate begins
- `qubit_count`: the number of qubits in the superposition this gate will act on
"""
function gate_extension(gate::Gate, qubit_index::Int64, qubit_count::Int64)
  qubit_index > qubit_count && throw(DomainError(qubit_index, "index cannot be greater than count"))
  matrix::Matrix{Complex{Float64}} = ((1.0 + 0.0im)*I)(2^(qubit_index - 1))
  matrix = kron(matrix, gate.matrix)
  identity = ((1.0 + 0.0im)*I)(2^(qubit_count - gate.superposed_qubits_required - qubit_index + 1))
  Gate(kron(matrix, identity))
end

"""
    gate_control(gate, control_qubits, qubit_index, qubit_count)

creates an n-controlled single qubit gate. acts on an arbitrary
qubit in a superposition and is controlled by n arbitrary qubits.

# Arguments
- `gate`: the underlying single qubit gate
- `control_qubits`: an array of integers selecting control qubits
- `qubit_index`: an integer selecting the qubit to which the single-qubit
gate is applied
- `qubit_count`: the number of qubits in the superposition this gate will
act on
"""
function gate_control(gate::Gate, control_qubits::Array{Int64}, qubit_index::Int64, qubit_count::Int64)
  qubit_index in control_qubits && throw(DomainError((control_qubit, qubit_index), "control and target qubits must differ"))
  
  outer_size = 2^qubit_count
  matrix::Matrix{Complex{Float64}} = Matrix{Complex{Float64}}(undef, (outer_size, outer_size))
  
  mask = 0
  for control_qubit in control_qubits
    mask += (1 << (qubit_count - control_qubit))
  end
  output_mask = (1 << (qubit_count - qubit_index))

  qubits::Matrix{Complex{Float64}} = Matrix{Complex{Float64}}(undef, (qubit_count, 2))
  for n in 1:outer_size
    for i in 1:qubit_count
      if i == qubit_index && ((n - 1) & mask) == mask
        qubits[i:i, :] = (((n - 1) & output_mask) == 0) ? gate.matrix[:, 1:1]' : gate.matrix[:, 2:2]'
      else
        qubits[i:i, :] = ((n - 1) & 2^(qubit_count - i)) == 0 ? gate_i.matrix[:, 1:1]' : gate_i.matrix[:, 2:2]'
      end
    end
    matrix[:, n:n] = qubit_tensor_product(qubits)
  end

  Gate(matrix)
end

"""
    gate_cx(control_qubit, qubit_index, qubit_count)

controlled pauli x gate.

# Arguments
- `control_qubit`: the index of the qubit controlling the pauli x operation
- `qubit_index`: the index of the qubit upon which to perform the pauli x operation
- `qubit_count`: the total number of qubits the resultant gate acts on
"""
function gate_cx(control_qubit::Int64, qubit_index::Int64, qubit_count::Int64)
  gate_control(gate_x, [control_qubit], qubit_index, qubit_count)
end

"""
    gate_cnx(control_qubits, qubit_index, qubit_count)

n-controlled pauli x gate.

# Arguments
- `control_qubits`: the indices of the qubits controlling the pauli x operation
- `qubit_index`: the index of the qubit upon which to perform the pauli x operation
- `qubit_count`: the total number of qubits the resultant gate acts on
"""
function gate_cnx(control_qubits::Array{Int64}, qubit_index::Int64, qubit_count::Int64)
  gate_control(gate_x, control_qubits, qubit_index, qubit_count)
end

# this special function creates a controlled-n gate efficiently, by utilizing the top n qubits
function gate_multi_control(gate::Gate, control_count::Int64, qubit_count::Int64)
  gate_qubit_count::Int64 = gate.superposed_qubits_required
  gate_qubit_count + control_count == qubit_count || throw(DomainError(qubit_count, "incorrect qubit count for gate + controls"))
  
  matrix::Matrix{Complex{Float64}} = ((1.0 + 0.0im)*I)(2^qubit_count)
  
  range = (2^qubit_count - 2^gate_qubit_count + 1):(2^qubit_count)
  matrix[range, range] = gate.matrix

  Gate(matrix)
end

"""
    gate_fourier_transform(qubit_count)

the quantum fourier transform.
https://en.wikipedia.org/wiki/Quantum_Fourier_transform

# Arguments
- `qubit_count`: the number of qubits the gate operates on
"""
function gate_fourier_transform(qubit_count::Int64)
  gate_size = 2^qubit_count
  ω = ℯ^(2im * pi / gate_size)

  matrix::Matrix{Complex{Float64}} = Matrix{Complex{Float64}}(undef, (gate_size, gate_size))
  for i in 1:gate_size
    for j in i:gate_size
      value = (ω^((i - 1) * (j - 1)))
      matrix[i, j] = value
      matrix[j, i] = value
    end
  end

  Gate(matrix / √gate_size)
end

"""
    gate_invert(gate)

inverts a unitary matrix quickly, taking advantage of the fact that the
conjugate transpose is the inverse of a unitary matrix

# Arguments
- `gate`: the gate to invert
"""
function gate_invert(gate::Gate)
  Gate(Array(gate.matrix'))
end

"""
    Measurement(bits_to_output, qubits_to_measure, sample_size)

creates a measurement component for incorporation in a `Circuit`.
`bits_to_output` and `qubits_to_measure` must be equal in length and
contain no duplicates. it is fine for there to be an intersection
of the to arrays, as long as each array has unique values.

# Arguments
- `bits_to_output`: the bits to write in the classical output register.
order must correspond to `qubits_to_measure` order.
- `qubits_to_measure`: the qubits to measure, in the order of measurement.
- `sample_size`: the number of samples to perform when measuring
"""
struct Measurement
  qubits_to_measure::Array{Int64, 1}
  bits_to_output::Array{Int64, 1}
  sample_size::Int64
  samples::Array{Int64, 1}

  function Measurement(bits_to_output::Array{Int64, 1}, qubits_to_measure::Array{Int64, 1}, sample_size::Int64 = 1024)
    length(qubits_to_measure) == length(bits_to_output) || throw(DomainError(length(bits_to_output), "tuples must have same number of unique elements"))
    qubits_to_measure = unique(qubits_to_measure)
    length(qubits_to_measure) == length(bits_to_output) || throw(DomainError(length(bits_to_output), "tuples must have same number of unique elements"))
    bits_to_output = unique(bits_to_output)
    length(qubits_to_measure) == length(bits_to_output) || throw(DomainError(length(bits_to_output), "tuples must have same number of unique elements"))

    new(qubits_to_measure, bits_to_output, sample_size, Array{Int64}(undef, sample_size))
  end
end

"""
    measure_superposition(superposition, classical_register, measurement)

i believe this function may not collapse things correctly, in that it may yield
a subset of the results that a physical quantum computer would.

# Arguments
- `superposition`: the superposition being measured
- `classical_register`: the output register
- `measurement`: details of the measurement
"""
function measure_superposition(superposition::Superposition, classical_register::ClassicalRegister, measurement::Measurement)
  measurement_qubit_count::Int64 = length(measurement.qubits_to_measure)
  register_qubit_count::Int64 = log2(length(superposition.state))
  probability_of_ones::Array{Float64} = zeros(Float64, measurement_qubit_count)

  for value in 0:(length(superposition.state) - 1)
    for i in 1:measurement_qubit_count
      qubit_to_measure = measurement.qubits_to_measure[i]
      exponent = register_qubit_count - qubit_to_measure
      if 2^exponent & value != 0
        probability_of_ones[i] += abs(superposition.state[value + 1] ^ 2)
      end
    end
  end

  for j in 1:measurement.sample_size
    for i in 1:measurement_qubit_count
      bit_to_output = measurement.bits_to_output[i]
      mask = 2^(classical_register.width - bit_to_output)
      if probability_of_ones[i] >= rand()
        classical_register.value |= mask
      else
        classical_register.value &= ~mask
      end
    end
    measurement.samples[j] = classical_register.value
  end

  counts = countmap(measurement.samples)
  max_value = first(counts)[1]
  max_count = first(counts)[2]
  for (value, count) in counts
    if count > max_count
      max_value = value
      max_count = count
    end
  end

  classical_register.value = max_value

  for i in 1:measurement_qubit_count
    bit_output = measurement.bits_to_output[i]
    output_mask::Int64 = 2^(classical_register.width - bit_output)
    bit_set = (classical_register.value & output_mask != 0)
    for value::Int64 in 0:(length(superposition.state) - 1)
      qubit_to_measure = measurement.qubits_to_measure[i]
      input_mask::Int64 = 2^(register_qubit_count - qubit_to_measure)
      if bit_set
        if value & input_mask == 0
          index = (value | input_mask) + 1
          if superposition.state[index] == 0
            # do we need to worry about phase here?
            superposition.state[index] = superposition.state[value + 1]
          else
            # the resulting phase ignores the component being zeroed out, is this correct?
            scale = √(abs(superposition.state[index])^2 + abs(superposition.state[value + 1])^2) / abs(superposition.state[index])
            superposition.state[index] *= scale
          end
          superposition.state[value + 1] = 0
        end
      else
        if value & input_mask != 0
          index = (value & ~input_mask) + 1
          if superposition.state[index] == 0
            # do we need to worry about phase here?
            superposition.state[index] = superposition.state[value + 1]
          else
            # the resulting phase ignores the component being zeroed out, is this correct?
            scale = √(abs(superposition.state[index])^2 + abs(superposition.state[value + 1])^2) / abs(superposition.state[index])
            superposition.state[index] *= scale
          end
          superposition.state[value + 1] = 0
        end
      end
    end
  end
end

"""
    HybridComponent(executor, arguments)

a component that can be added to a `Circuit` that is capable of
controlling quantum gates and subcircuits with classical logic.

# Arguments
- `executor`: a function that accepts circuit application arguments (see
below for signature) and runs code that applies quantum `Gates` and `Circuits`
intelligently
- `arguments`: fixed arguments known at the time of `Circuit` composition
"""
struct HybridComponent
  execute
  arguments::Array

  # executor should be a function with the signature
  # (superposition::QuantumComputer.Superposition, classical_register::QuantumComputer.ClassicalRegister, arguments...)
  function HybridComponent(executor, arguments)
    new(executor, arguments)
  end
end

"""
    Circuit()

a quantum circuit. technically this is just an ordered collection of
components.
"""
struct Circuit
  components::Array{Union{Circuit, Gate, Measurement}, 1}

  function Circuit()
    new(Array{Union{Circuit, Gate, Measurement}}(undef, 0))
  end
end

"""
    add_gate_to_circuit!(circuit, gate)

adds a `Gate` to a `Circuit`.

# Arguments
- `circuit`: the `Circuit`
- `gate`: the `Gate`
"""
function add_gate_to_circuit!(circuit::Circuit, gate::Gate)
  push!(circuit.components, gate)
end

"""
    add_measurement_to_circuit!(circuit, measurement)

adds a `Measurement` to a `Circuit`.

# Arguments
- `circuit`: the `Circuit`
- `measurement`: the `Measurement`
"""
function add_measurement_to_circuit!(circuit::Circuit, measurement::Measurement)
  push!(circuit.components, measurement)
end

"""
    add_hybrid_component_to_circuit!(circuit, component)

adds a `HybridComponent` to a `Circuit`.

# Arguments
- `circuit`: the `Circuit`
- `component`: the `HybridComponent`
"""
function add_hybrid_component_to_circuit!(circuit::Circuit, component::HybridComponent)
  push!(circuit.components, component)
end

"""
    add_subcircuit_to_circuit!(circuit, subcircuit)

adds a `Circuit` to a another `Circuit` as a component.

# Arguments
- `circuit`: the parent `Circuit`
- `subcircuit`: the child `Circuit`
"""
function add_subcircuit_to_circuit!(circuit::Circuit, subcircuit::Circuit)
  push!(circuit.components, subcircuit)
end

"""
    apply_circuit_to_superposition!(superposition, circuit, classical_register)

applies the gates and measurements in a circuit to the superposition provided.
a classical register is required for measurement outputs.

# Arguments
- `superposition`: the `Superposition` to operate on
- `circuit`: the `Circuit` to apply
- `classical_register`: a `ClassicalRegister` capable of holding all circuit
measurement output bits
"""
function apply_circuit_to_superposition!(superposition::Superposition, circuit::Circuit, classical_register::ClassicalRegister = undef)
  for component in circuit.components
    if typeof(component) == Circuit
      # if we are recursing too deep our circuit is likely pretty poor
      # so this is reasonable
      apply_circuit_to_superposition!(superposition, component, classical_register)
    elseif typeof(component) == Measurement
      measure_superposition(superposition, classical_register, component)
    elseif typeof(component) == HybridComponent
      component.execute(superposition, classical_register, component.arguments...)
    else
      superposition.state = component.matrix * superposition.state
    end
  end
end

"""
    QuantumComputer.Circuits

a collection of simple quantum circuits.

# Example
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
"""
module Circuits

using LinearAlgebra
using ..QuantumComputer

"""
    constant_adder(n, qubit_count)

a quantum circuit that adds n to the superposition's value

# Arguments:
- `n`: the constant to add
- `qubit_count`: the number of qubits in the superposition (the adder will
also return the result modulo 2^qubit_count)
"""
function constant_adder(n::Int64, qubit_count::Int64)
  0 <= n < 2^qubit_count || throw(DomainError(n, ""))

  circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()
  qft::QuantumComputer.Gate = QuantumComputer.gate_fourier_transform(qubit_count)
  QuantumComputer.add_gate_to_circuit!(circuit, qft)
  for i::Int64 in qubit_count:-1:1
    for k::Int64 in 1:i
      if (1 << (i - k)) & n != 0
        gate::QuantumComputer.Gate = QuantumComputer.gate_extension(QuantumComputer.gate_p(pi / 2^(k - 1)), i, qubit_count)
        QuantumComputer.add_gate_to_circuit!(circuit, gate)
      end
    end
  end
  inverse_qft::QuantumComputer.Gate = QuantumComputer.gate_invert(qft)
  QuantumComputer.add_gate_to_circuit!(circuit, inverse_qft)

  circuit
end

"""
    period_finding_for_11x_mod_15()

a 5 qubit circuit that returns a value related to the period of
f(x) = 11^x mod 15. this effectively allows implementation of Shor's
algorithm for this specific case.
"""
function period_finding_for_11x_mod_15()
  qubit_count = 5

  circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()
  h_matrix = QuantumComputer.gate_h.matrix
  gate_h::QuantumComputer.Gate = QuantumComputer.Gate(kron(((1.0+0.0im)*I)(4), kron(h_matrix, kron(h_matrix, h_matrix))))
  QuantumComputer.add_gate_to_circuit!(circuit, gate_h)
  gate_cx::QuantumComputer.Gate = QuantumComputer.gate_cx(3, 4, qubit_count)
  QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
  gate_cx = QuantumComputer.gate_cx(3, 5, qubit_count)
  QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
  gate_h = QuantumComputer.gate_extension(QuantumComputer.gate_h, 4, qubit_count)
  QuantumComputer.add_gate_to_circuit!(circuit, gate_h)
  gate_cp = QuantumComputer.gate_control(QuantumComputer.gate_p(pi/2), [4], 5, qubit_count)
  QuantumComputer.add_gate_to_circuit!(circuit, gate_cp)
  gate_h = QuantumComputer.gate_extension(QuantumComputer.gate_h, 5, qubit_count)
  QuantumComputer.add_gate_to_circuit!(circuit, gate_h)

  # uncompute ancilla
  gate_cx = QuantumComputer.gate_cx(3, 2, qubit_count)
  QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
  gate_cx = QuantumComputer.gate_cx(3, 1, qubit_count)
  QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
  # done uncomputing

  gate_cp = QuantumComputer.gate_control(QuantumComputer.gate_p(pi/4), [4], 3, qubit_count)
  QuantumComputer.add_gate_to_circuit!(circuit, gate_cp)
  gate_cp = QuantumComputer.gate_control(QuantumComputer.gate_p(pi/2), [5], 3, qubit_count)
  QuantumComputer.add_gate_to_circuit!(circuit, gate_cp)

  circuit
end

end # module Circuits

end # module QuantumComputer
