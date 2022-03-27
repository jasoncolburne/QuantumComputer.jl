module QuantumComputer

using LinearAlgebra

function qubit_tensor_product_thread(next_product::Array{Complex{Float64}, 1}, next_basis::Array{Complex{Float64}, 1}, product::Array{Complex{Float64}, 1}, product_size::Int64, thread_number::Int64, thread_count::Int64)
  width::Int64, padded_count::Int64 = divrem(product_size, thread_count)

  if padded_count < thread_number
    initial_value = (width + 1) * padded_count + width * (thread_number - padded_count - 1) + 1
  else
    width += 1
    initial_value = width * (thread_number - 1) + 1
  end

  for i::Int64 in initial_value:min(initial_value + width - 1, product_size)
    next_product[2i - 1] = next_basis[1] * product[i]
    next_product[2i] = next_basis[2] * product[i]
  end
end

# this custom function computes a tensor product on the rows of
# a matrix of qubits (all are considered orthogonal)
function qubit_tensor_product(qubits::Matrix{Complex{Float64}})::Array{Complex{Float64}, 1}
  qubit_count::Int64 = size(qubits, 1)
  thread_count::Int64 = Threads.nthreads()

  product::Array{Complex{Float64}, 1} = qubit_count > 1 ? qubits[1, :] : throw(DomainError(qubits, "at least 2 qubits required"))
  for i::Int64 in 2:qubit_count
    next_basis::Array{Complex{Float64}, 1} = qubits[i, :]
    product_size::Int64 = length(product)
    next_product::Array{Complex{Float64}, 1} = Array{Complex{Float64}}(undef, 2 * product_size)
    
    Threads.@threads for thread_number::Int64 in 1:min(thread_count, product_size)
      qubit_tensor_product_thread(next_product, next_basis, product, product_size, thread_number, thread_count)
    end
    
    product = next_product
  end

  product
end

mutable struct ClassicalRegister
  width::Int64
  value::Int64

  function ClassicalRegister(width::Int64, value::Int64 = 0)
    value < 2^width || throw(DomainError(value, "value does not fit in register"))

    new(width, value)
  end
end

mutable struct Register
  qubits::Matrix{Complex{Float64}}

  function Register(qubit_count::Int64, value::Int64 = 0)
    qubit_count > 32 && throw(DomainError(qubit_count, "you are truly mad"))
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

mutable struct Superposition
  state::Array{Complex{Float64}, 1}

  function Superposition(qubits::Matrix{Complex{Float64}})
    new(qubit_tensor_product(qubits))
  end

  function Superposition(state::Array{Complex{Float64}, 1})
    new(state)
  end
end

function qubits_operated_on_by_unitary_matrix(matrix::Matrix)::Int64
  dimensions = size(matrix)
  dimensions[1] == dimensions[2] || throw(DomainError(matrix, "matrix must be square"))
  log2(dimensions[1])
end

struct Gate
  matrix::Matrix{Complex{Float64}}
  superposed_qubits_required::Int64

  function Gate(matrix::Matrix)
    new(matrix, qubits_operated_on_by_unitary_matrix(matrix))
  end
end

# single qubit gates
# https://qiskit.org/textbook/ch-states/single-qubit-gates.html
####################

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
gate_root_swap = Gate([1 0 0 0; 0 (1.0+im)/2 (1.0-im)/2 0; 0 (1.0-im)/2 (1.0+im)/2 0; 0 0 0 1])

# https://quantumcomputing.stackexchange.com/questions/9181/swap-gate-on-2-qubits-in-3-entangled-qubit-system
function ket_bra(i::Int64, j::Int64)::Matrix{Complex{Float64}}
  (i > 1 || i < 0 || j > 1 || j < 0) && throw(DomainError((i, j), "arguments must be in range [0, 1]"))

  matrix::Matrix{Complex{Float64}} = zeros(Complex{Float64}, (2, 2))
  matrix[[(j << 1) + i + 1]] = [1.0 + 0.0im]
  matrix
end

function gate_swap(qubit_a_index::Int64, qubit_b_index::Int64, qubit_count::Int64)
  identity = ((1.0 + 0.0im)*I)(2)

  matrix_size = 2^qubit_count
  matrix::Matrix{Complex{Float64}} = zeros(Complex{Float64}, (matrix_size, matrix_size))
  matrices::Array{Matrix{Complex{Float64}}} = Array{Matrix{Complex{Float64}}}(undef, 0)
  
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

function gate_extension(gate::Gate, qubit_index::Int64, qubit_count::Int64)
  qubit_index > qubit_count && throw(DomainError(qubit_index, "index cannot be greater than count"))
  matrix::Matrix{Complex{Float64}} = ((1.0 + 0.0im)*I)(2^(qubit_index - 1))
  matrix = kron(matrix, gate.matrix)
  identity = ((1.0 + 0.0im)*I)(2^(qubit_count - gate.superposed_qubits_required - qubit_index + 1))
  Gate(kron(matrix, identity))
end

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

function gate_cx(control_qubit::Int64, qubit_index::Int64, qubit_count::Int64)
  gate_control(gate_x, [control_qubit], qubit_index, qubit_count)
end

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

# way faster for a unitary matrix than inv()
function gate_invert(gate::Gate)
  Gate(Array(gate.matrix'))
end

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

      mask::Int64 = 2^(classical_register.width - bit_to_output)
      if probability_of_ones[i] >= rand()
        classical_register.value |= mask
      else
        classical_register.value &= ~mask
      end
    end
    measurement.samples[j] = classical_register.value
  end
end

struct Circuit
  components::Array{Union{Circuit, Gate, Measurement}, 1}

  function Circuit()
    new(Array{Union{Circuit, Gate, Measurement}}(undef, 0))
  end
end

function add_gate_to_circuit!(circuit::Circuit, gate::Gate)
  push!(circuit.components, gate)
end

function add_measurement_to_circuit!(circuit::Circuit, measurement::Measurement)
  push!(circuit.components, measurement)
end

function add_subcircuit_to_circuit!(circuit::Circuit, subcircuit::Circuit)
  push!(circuit.components, subcircuit)
end

function apply_circuit_to_superposition!(superposition::Superposition, circuit::Circuit, classical_register::ClassicalRegister = undef)
  for component in circuit.components
    if typeof(component) == Circuit
      # if we are recursing too deep our circuit is likely pretty poor
      # so this is reasonable
      apply_circuit_to_superposition!(superposition, component, classical_register)
    elseif typeof(component) == Measurement
      measure_superposition(superposition, classical_register, component)
    else
      superposition.state = component.matrix * superposition.state
    end
  end
end

module Circuits

using LinearAlgebra
using ..QuantumComputer

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

function period_finding(a::Int64, n::Int64)
end

end

end # module
