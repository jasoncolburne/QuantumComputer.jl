"""
    using QuantumComputer

a quantum computer simulator. built out of boredom and intrigue. if you are new to julia, be warned - array indices start at 1. throughout this module, qubit index 1 is considered the most significant qubit. in typical circuit diagrams, it represents the top line.

i am fully aware that i 'over-typed' this module. i am new to julia and was curious about performance differences but i assume the compiler will sort everything out in many of the cases i am being explicit.

what this module can do:
- create registers, superpositions, gates and circuits
- apply circuits to superpositions
- add gates to circuits
- add measurements to circuits
- add circuits to other circuits as subcircuits
- add `HybridComponents` to a circuit, which are components where there is a classical interaction with the quantum portion of the simulator
- measure arbitrary qubits from a superposition
- generate a variety of single- and multi-qubit static and dynamic gates
- generate some circuits
- construct a gate from a circuit
- cache gates on disk

what this module cannot do:
- decompose a gate into a circuit of smaller gates (planned, long term)
- print a graphical representation of the circuit (planned, relatively simple but also tedious so not planned for near future)
- use multiple registers in one circuit (right now you get one quantum and one classical)
- output bloch sphere animations of algorithms running (planned)

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
using Serialization
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

        for i::Int64 = 1:qubit_count
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

a superposition is a representation of the probabilities of all possible states of a vector of qubits, typically a register (but there is no hard requirement for a register as input in this simulator)
"""
mutable struct Superposition
    state::Array{Complex{Float64},1}

    function Superposition(qubits::Matrix{Complex{Float64}})
        new(qubit_tensor_product(qubits))
    end

    function Superposition(state::Array{Complex{Float64},1})
        new(state)
    end
end

"""
    qubit_tensor_product(qubits)

this custom function computes a tensor product on the n rows of a matrix of 2-column qubits. this function is threaded, and the core of the algorithm follows in the next function. this works by iteratively applying the tensor product qubit by qubit.

# Arguments
- `qubits`: a matrix of qubits, in the same arrangement as a `Register`
"""
function qubit_tensor_product(qubits::Matrix{Complex{Float64}})::Array{Complex{Float64},1}
    qubit_count::Int64 = size(qubits, 1)
    thread_count::Int64 = Threads.nthreads()

    product::Array{Complex{Float64},1} =
        qubit_count > 1 ? qubits[1, :] :
        throw(DomainError(qubits, "at least 2 qubits required"))
    for i::Int64 = 2:qubit_count
        basis::Array{Complex{Float64},1} = qubits[i, :]
        product_size::Int64 = length(product)
        next_product::Array{Complex{Float64},1} =
            Array{Complex{Float64}}(undef, 2 * product_size)

        # coming back to this code, i was confused why we're iterating to 
        # min(thread_count, product_size) but passing the child function thread_count
        # turns out that since we need to do some computations with product size in
        # the child function, we end up handling this there.
        Threads.@threads for thread_number::Int64 = 1:min(thread_count, product_size)
            qubit_tensor_product_thread(
                next_product,
                basis,
                product,
                product_size,
                thread_number,
                thread_count,
            )
        end

        product = next_product
    end

    product
end

"""
    qubit_tensor_product_thread(next_product, basis, product, product_size, thread_number, thread_count)

the threaded portion from our qubit tensor product algorithm. this will compute a slice of the current tensor product expansion from the parent algorithm.
this fills in a portion of `next_product` by computing a partial tensor product of `basis` and `product`.

# Arguments
- `next_product`: the product (state) that is currently being computed
- `basis`: the qubit being applied to the previous product
- `product`: the previous product
- `product_size`: the number of states represented by `product`
- `thread_number`: the algorithm-assigned thread number of this thread
- `thread_count`: the total number of threads
"""
function qubit_tensor_product_thread(
    next_product::Array{Complex{Float64},1},
    basis::Array{Complex{Float64},1},
    product::Array{Complex{Float64},1},
    product_size::Int64,
    thread_number::Int64,
    thread_count::Int64,
)
    width::Int64, padded_count::Int64 = divrem(product_size, thread_count)

    if padded_count < thread_number
        initial_value =
            (width + 1) * padded_count + width * (thread_number - padded_count - 1) + 1
    else
        width += 1
        initial_value = width * (thread_number - 1) + 1
    end

    for i::Int64 = initial_value:min(initial_value + width - 1, product_size)
        next_product[2i-1] = basis[1] * product[i]
        next_product[2i] = basis[2] * product[i]
    end
end

# gates
#######

"""
    Gate(matrix)

a gate operating on a state of size 2^n is really just a 2^n x 2^n unitary matrix. this can be added to a `Circuit`.

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
gate_h = Gate([1 1; 1 -1] / ???2)
# phase
gate_p(??::Float64) = Gate([1 0; 0 ???^(??*im)])
# identity
gate_i = Gate([1 0; 0 1])
# pi / 2 and pi / 4 phase rotations
gate_s = gate_p(pi / 2)
gate_sdg = gate_p(-pi / 2)
gate_t = gate_p(pi / 4)
gate_tdg = gate_p(-pi / 4)
# u, the most general gate
gate_u(??::Float64, ??::Float64, ??::Float64) =
    Gate([cos(?? / 2) -sin(?? / 2)*???^(??*im); sin(?? / 2)*???^(??*im) cos(?? / 2)*???^((??+??)*im)])

# multi-qubit gates
###################

"""
[wikipedia](https://en.wikipedia.org/wiki/Quantum_logic_gate#Square_root_of_swap_gate)
"""
gate_root_swap =
    Gate([1 0 0 0; 0 (1.0+im)/2 (1.0-im)/2 0; 0 (1.0-im)/2 (1.0+im)/2 0; 0 0 0 1])

"""
    gate_swap(qubit_a_index, qubit_b_index, qubit_count)

a gate that swaps two arbitrary qubits in a superposition [explanation](https://quantumcomputing.stackexchange.com/questions/9181/swap-gate-on-2-qubits-in-3-entangled-qubit-system)

# Arguments
- `qubit_a_index`: the index of the first qubit to swap
- `qubit_b_index`: the index of the second qubit to swap
- `qubit_count`: the total number of qubits in the superposition
this gate will operate on
"""
function gate_swap(qubit_a_index::Int64, qubit_b_index::Int64, qubit_count::Int64)
    identity = ((1.0 + 0.0im) * I)(2)

    matrix_size = 2^qubit_count
    matrix::Matrix{Complex{Float64}} = zeros(Complex{Float64}, (matrix_size, matrix_size))
    matrices::Array{Matrix{Complex{Float64}},1} = Array{Matrix{Complex{Float64}}}(undef, 0)

    for i = 1:qubit_count
        push!(matrices, identity)
    end

    for i = 0:1
        for j = 0:1
            matrices[qubit_a_index] = ket_bra(i, j)
            matrices[qubit_b_index] = ket_bra(j, i)
            partial_sum = matrices[1]
            for k = 2:qubit_count
                partial_sum = kron(partial_sum, matrices[k])
            end
            matrix += partial_sum
        end
    end

    Gate(matrix)
end

"""
    ket_bra(i, j)

check the stackexchange post on the previous function to understand what is happening here.

# Arguments
- `i`: first value
- `j`: second value
"""
function ket_bra(i::Int64, j::Int64)::Matrix{Complex{Float64}}
    (i > 1 || i < 0 || j > 1 || j < 0) &&
        throw(DomainError((i, j), "arguments must be in range [0, 1]"))

    matrix::Matrix{Complex{Float64}} = zeros(Complex{Float64}, (2, 2))
    matrix[[(j << 1) + i + 1]] = [1.0 + 0.0im]
    matrix
end

function gate_expansion(gate::Gate, qubit_count::Int64)
    matrix_single = gate.matrix
    matrix = matrix_single
    for i in 2:qubit_count
         matrix = kron(matrix, matrix_single)
    end
    QuantumComputer.Gate(matrix)
end

"""
    gate_extension(gate, qubit_index, qubit_count)

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
    qubit_index > qubit_count &&
        throw(DomainError(qubit_index, "index cannot be greater than count"))
    matrix::Matrix{Complex{Float64}} = ((1.0 + 0.0im) * I)(2^(qubit_index - 1))
    matrix = kron(matrix, gate.matrix)
    identity = ((1.0 + 0.0im) * I)(
        2^(qubit_count - gate.superposed_qubits_required - qubit_index + 1),
    )
    Gate(kron(matrix, identity))
end

"""
    gate_control(gate, control_qubits, qubit_index, qubit_count)

creates an n-controlled single qubit gate. acts on an arbitrary qubit in a superposition and is controlled by n arbitrary qubits.

# Arguments
- `gate`: the underlying single qubit gate
- `control_qubits`: an array of integers selecting control qubits
- `qubit_index`: an integer selecting the qubit to which the single-qubit gate is applied
- `qubit_count`: the number of qubits in the superposition this gate will act on
"""
function gate_control(
    gate::Gate,
    control_qubits::Array{Int64},
    qubit_index::Int64,
    qubit_count::Int64,
)
    qubit_index in control_qubits && throw(
        DomainError((control_qubits, qubit_index), "control and target qubits must differ"),
    )

    outer_size = 2^qubit_count
    matrix::Matrix{Complex{Float64}} =
        Matrix{Complex{Float64}}(undef, (outer_size, outer_size))

    mask = 0
    for control_qubit in control_qubits
        mask += (1 << (qubit_count - control_qubit))
    end
    output_mask = (1 << (qubit_count - qubit_index))

    qubits::Matrix{Complex{Float64}} = Matrix{Complex{Float64}}(undef, (qubit_count, 2))
    for n = 1:outer_size
        for i = 1:qubit_count
            if i == qubit_index && ((n - 1) & mask) == mask
                qubits[i:i, :] =
                    (((n - 1) & output_mask) == 0) ? gate.matrix[:, 1:1]' :
                    gate.matrix[:, 2:2]'
            else
                qubits[i:i, :] =
                    ((n - 1) & 2^(qubit_count - i)) == 0 ? gate_i.matrix[:, 1:1]' :
                    gate_i.matrix[:, 2:2]'
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

"""
    gate_multi_control(gate, control_count, qubit_count)

this special function creates a controlled-n gate efficiently, by utilizing the top n qubits

# Arguments
- `gate`
"""
function gate_multi_control(gate::Gate, control_count::Int64, qubit_count::Int64)
    gate_qubit_count::Int64 = gate.superposed_qubits_required
    gate_qubit_count + control_count == qubit_count ||
        throw(DomainError(qubit_count, "incorrect qubit count for gate + controls"))

    matrix::Matrix{Complex{Float64}} = ((1.0 + 0.0im) * I)(2^qubit_count)

    range = (2^qubit_count-2^gate_qubit_count+1):(2^qubit_count)
    matrix[range, range] = gate.matrix

    Gate(matrix)
end

"""
    gate_fourier_transform(qubit_count)

the quantum fourier transform. [wikipedia](https://en.wikipedia.org/wiki/Quantum_Fourier_transform)

# Arguments
- `qubit_count`: the number of qubits the gate operates on
"""
function gate_fourier_transform(qubit_count::Int64)
    gate_size = 2^qubit_count
    ?? = ???^(2im * pi / gate_size)

    matrix::Matrix{Complex{Float64}} =
        Matrix{Complex{Float64}}(undef, (gate_size, gate_size))
    for i = 1:gate_size
        for j = i:gate_size
            value = (??^((i - 1) * (j - 1)))
            matrix[i, j] = value
            matrix[j, i] = value
        end
    end

    Gate(matrix / ???gate_size)
end

"""
    gate_invert(gate)

inverts a unitary matrix quickly, taking advantage of the fact that the conjugate transpose is the inverse of a unitary matrix

# Arguments
- `gate`: the gate to invert
"""
function gate_invert(gate::Gate)
    Gate(Array(gate.matrix'))
end

"""
    Measurement(bits_to_output, qubits_to_measure, sample_size)

creates a measurement component for incorporation in a `Circuit`. `bits_to_output` and `qubits_to_measure` must be equal in length and contain no duplicates. it is fine for there to be an intersection of the to arrays, as long as each array has unique values.

# Arguments
- `bits_to_output`: the bits to write in the classical output register. order must correspond to `qubits_to_measure` order.
- `qubits_to_measure`: the qubits to measure, in the order of measurement.
- `sample_size`: the number of samples to perform when measuring
"""
struct Measurement
    qubits_to_measure::Array{Int64,1}
    bits_to_output::Array{Int64,1}
    sample_size::Int64
    samples::Array{Int64,1}

    function Measurement(
        bits_to_output::Array{Int64,1},
        qubits_to_measure::Array{Int64,1},
        sample_size::Int64 = 1024,
    )
        length(qubits_to_measure) == length(bits_to_output) || throw(
            DomainError(
                length(bits_to_output),
                "tuples must have same number of unique elements",
            ),
        )
        qubits_to_measure = unique(qubits_to_measure)
        length(qubits_to_measure) == length(bits_to_output) || throw(
            DomainError(
                length(bits_to_output),
                "tuples must have same number of unique elements",
            ),
        )
        bits_to_output = unique(bits_to_output)
        length(qubits_to_measure) == length(bits_to_output) || throw(
            DomainError(
                length(bits_to_output),
                "tuples must have same number of unique elements",
            ),
        )

        new(
            qubits_to_measure,
            bits_to_output,
            sample_size,
            Array{Int64,1}(undef, sample_size),
        )
    end
end

# using Printf

"""
    measure_superposition(superposition, classical_register, measurement)

measure some qubits and store the result in a classical register

# Arguments
- `superposition`: the superposition being measured
- `classical_register`: the output register
- `measurement`: details of the measurement
"""
function measure_superposition(
    superposition::Superposition,
    classical_register::ClassicalRegister,
    measurement::Measurement,
)
    measurement_qubit_count::Int64 = length(measurement.qubits_to_measure)
    register_qubit_count::Int64 = log2(length(superposition.state))
    probability_of_ones::Array{Float64} = zeros(Float64, measurement_qubit_count)

    for value = 0:(length(superposition.state)-1)
        for i = 1:measurement_qubit_count
            qubit_to_measure = measurement.qubits_to_measure[i]
            exponent = register_qubit_count - qubit_to_measure
            if 2^exponent & value != 0
                probability_of_ones[i] += abs(superposition.state[value+1])^2
            end
        end
    end

    initial_value = classical_register.value
    for j = 1:measurement.sample_size
        measurement.samples[j] = initial_value
        for i = 1:measurement_qubit_count
            bit_to_output = measurement.bits_to_output[i]
            mask = 2^(classical_register.width - bit_to_output)
            # should this be >= or >? it surely doesn't matter but i'd like to be correct
            if probability_of_ones[i] >= rand()
                measurement.samples[j] |= mask
            else
                measurement.samples[j] &= ~mask
            end
        end
    end

    counts = countmap(measurement.samples)
    most_common_value = first(counts)[1]
    max_count = first(counts)[2]
    for (value, count) in counts
        if count > max_count
            most_common_value = value
            max_count = count
        end
    end

    classical_register.value = most_common_value
    # @printf(
    #     "classical register value is now %s\n",
    #     string(classical_register.value, base = 2, pad = classical_register.width)
    # )

    for i = 1:measurement_qubit_count
        bit_output = measurement.bits_to_output[i]
        output_mask::Int64 = 2^(classical_register.width - bit_output)
        bit_set = (classical_register.value & output_mask != 0)
        for value::Int64 = 0:(length(superposition.state)-1)
            qubit_to_measure = measurement.qubits_to_measure[i]
            input_mask::Int64 = 2^(register_qubit_count - qubit_to_measure)
            if bit_set
                if value & input_mask == 0
                    index = (value | input_mask) + 1
                    if superposition.state[index] == 0
                        # do we need to worry about phase here?
                        superposition.state[index] = superposition.state[value+1]
                    else
                        # the resulting phase ignores the component being zeroed out, is this correct?
                        scale =
                            ???(
                                abs(superposition.state[index])^2 +
                                abs(superposition.state[value+1])^2,
                            ) / abs(superposition.state[index])
                        superposition.state[index] *= scale
                    end
                    superposition.state[value+1] = 0
                end
            else
                if value & input_mask != 0
                    index = (value & ~input_mask) + 1
                    if superposition.state[index] == 0
                        # do we need to worry about phase here?
                        superposition.state[index] = superposition.state[value+1]
                    else
                        # the resulting phase ignores the component being zeroed out, is this correct?
                        scale =
                            ???(
                                abs(superposition.state[index])^2 +
                                abs(superposition.state[value+1])^2,
                            ) / abs(superposition.state[index])
                        superposition.state[index] *= scale
                    end
                    superposition.state[value+1] = 0
                end
            end
        end
    end
end

"""
    HybridComponent(executor, arguments)

a component that can be added to a `Circuit` that is capable of controlling quantum gates and subcircuits with classical logic.

# Arguments
- `executor`: a function that accepts circuit application arguments (see below for signature) and runs code that applies quantum `Gates` and `Circuits` intelligently
- `arguments`: fixed arguments known at the time of `Circuit` composition
"""
struct HybridComponent
    execute::Any
    arguments::Array

    # executor should be a function with the signature
    # (superposition::QuantumComputer.Superposition, classical_register::QuantumComputer.ClassicalRegister, arguments...)
    function HybridComponent(executor, arguments)
        new(executor, arguments)
    end
end

"""
    Circuit()

a quantum circuit. technically this is just an ordered collection of components.
"""
struct Circuit
    components::Array{Union{Circuit,Gate,HybridComponent,Measurement},1}

    function Circuit()
        new(Array{Union{Circuit,Gate,HybridComponent,Measurement}}(undef, 0))
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

cache_base = nothing
cached_gates = Dict{String,Gate}()

function gate_set_cache_base(base)
    global cache_base = [base]
end

function gate_load_from_cache(cache_path)
    key = join(cache_path, "/")
    key in keys(cached_gates) && return cached_gates[key]
    cache_base == nothing && return nothing

    # after this pop, we need to make it to the push at the end of the function
    # so don't return early for some reason
    filename = pop!(cache_path)

    path = join(vcat(cache_base, cache_path), "/")
    run(`mkdir -p $path`)
    cd(path)

    gate = nothing
    try
        gate = deserialize(string(filename, ".qg"))
        cached_gates[key] = gate
    catch
    end

    ancestor_path = [".." for _ in vcat(cache_base, cache_path)]
    path = join(ancestor_path, "/")
    cd(path)

    push!(cache_path, filename)

    return gate
end

function gate_save_to_cache(cache_path, gate)
    if cache_base != nothing
        filename = pop!(cache_path)

        path = join(vcat(cache_base, cache_path), "/")
        run(`mkdir -p $path`)
        cd(path)

        serialize(string(filename, ".qg"), gate)

        ancestor_path = [".." for _ in vcat(cache_base, cache_path)]
        path = join(ancestor_path, "/")
        cd(path)

        push!(cache_path, filename)
    end

    key = join(cache_path, "/")
    cached_gates[key] = gate
end

function gate_remove_from_cache(cache_path)
    if cache_base != nothing
        filename = pop!(cache_path)

        path = join(vcat(cache_base, cache_path), "/")
        run(`mkdir -p $path`)
        try
            cd(path)
            rm(string(filename, ".qg"))
        catch
        finally
            # this is pretty fragile but the mkdir above kind of saves us from edge cases
            ancestor_path = [".." for _ in vcat(cache_base, cache_path)]
            path = join(ancestor_path, "/")
            cd(path)
        end

        push!(cache_path, filename)
    end

    key = join(cache_path, "/")
    delete!(cached_gates, key)
end

function circuit_convert_to_gate(circuit::Circuit, cache_path::Array{String,1} = Array{String,1}(undef, 0))
    if length(cache_path) > 0
        gate = gate_load_from_cache(cache_path)
        if typeof(gate) != Nothing
            return gate
        end
    end

    gates = filter(component -> typeof(component) == Gate, circuit.components)
    rank = size(gates[1].matrix, 1)

    # initialize with identity
    matrix = ((1.0 + 0.0im) * I)(rank)

    for component in circuit.components
        if typeof(component) == Circuit
            # if we are recursing too deep our circuit is likely pretty poor
            # so this is reasonable
            matrix = circuit_convert_to_gate(component).matrix * matrix
        elseif typeof(component) == Measurement
            throw(
                DomainError(circuit, "cannot contain measurements when converting to gate"),
            )
        elseif typeof(component) == HybridComponent
            throw(
                DomainError(
                    circuit,
                    "cannot contain hybrid components when converting to gate",
                ),
            )
        else
            matrix = component.matrix * matrix
        end
    end

    gate = Gate(matrix)

    if length(cache_path) > 0
      gate_save_to_cache(cache_path, gate)
    end

    return gate
end

"""
    apply_circuit_to_superposition!(superposition, circuit, classical_register)

applies the gates and measurements in a circuit to the superposition provided. a classical register is required for measurement outputs.

# Arguments
- `superposition`: the `Superposition` to operate on
- `circuit`: the `Circuit` to apply
- `classical_register`: a `ClassicalRegister` capable of holding all circuit measurement output bits
"""
function apply_circuit_to_superposition!(
    superposition::Superposition,
    circuit::Circuit,
    classical_register::ClassicalRegister = undef,
)
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

a quantum circuit that adds `n` to the superposition's value (`mod 2^qubit_count`). based on Draper's circuit.

# Arguments:
- `n`: the constant to add
- `qubit_count`: the number of qubits in the superposition
"""
function constant_adder(n::Int64, qubit_count::Int64)
    0 <= n < 2^qubit_count || throw(DomainError(n, "must be less than 2^qubit_count"))

    circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()
    qft::QuantumComputer.Gate = QuantumComputer.gate_fourier_transform(qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, qft)
    subcircuit::QuantumComputer.Circuit = constant_adder_core(n, qubit_count)
    QuantumComputer.add_subcircuit_to_circuit!(circuit, subcircuit)
    inverse_qft::QuantumComputer.Gate = QuantumComputer.gate_invert(qft)
    QuantumComputer.add_gate_to_circuit!(circuit, inverse_qft)

    circuit
end

"""
    constant_adder_core(n, qubit_count)

a quantum circuit that adds `n` to the superposition's value (`mod 2^qubit_count`) without applying pre and post fourier transforms

# Arguments:
- `n`: the constant to add
- `qubit_count`: the number of qubits in the superposition
"""
function constant_adder_core(n::Int64, qubit_count::Int64)
    0 <= n < 2^qubit_count || throw(DomainError(n, "must be less than 2^qubit_count"))

    circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()
    for i::Int64 = qubit_count:-1:1
        for k::Int64 = 1:i
            if (1 << (i - k)) & n != 0
                gate::QuantumComputer.Gate = QuantumComputer.gate_extension(
                    QuantumComputer.gate_p(pi / 2^(k - 1)),
                    i,
                    qubit_count,
                )
                QuantumComputer.add_gate_to_circuit!(circuit, gate)
            end
        end
    end

    circuit
end

"""
    shor2n3_controlled_controlled_modular_adder(n, a, qubit_count)

a quantum circuit that adds `a` to the superposition's value (`mod n`). see the modular adder circuit in [this paper](https://arxiv.org/pdf/quant-ph/0205095.pdf).

# Arguments:
- `n`: the modulus
- `a`: the constant to add
- `no_cache`: if true, caching will not be employed when constructing gates
- `rebuild`: if true, this will destroy the relevant cached gates before regenerating them
"""
function shor2n3_controlled_controlled_modular_adder(n::Int64, a::Int64, no_cache = true, rebuild = false)
    n_qubit_count::Int64 = ceil(log2(n))
    qubit_count::Int64 = 2 * n_qubit_count + 3

    cache_path = no_cache ? Array{String,1}(undef, 0) : ["constant_adder_core", string(n_qubit_count + 1), string(a)]
    rebuild && !no_cache && gate_remove_from_cache(cache_path)
    gate_add_a = no_cache ? nothing : QuantumComputer.gate_load_from_cache(cache_path)
    if typeof(gate_add_a) == Nothing
        circuit_add_a::QuantumComputer.Circuit = constant_adder_core(a, n_qubit_count + 1)
        gate_add_a = QuantumComputer.circuit_convert_to_gate(circuit_add_a, cache_path)
    end
    gate_subtract_a::QuantumComputer.Gate = QuantumComputer.gate_invert(gate_add_a)

    offset::Int64 = qubit_count - gate_add_a.superposed_qubits_required - 1
    add_a::QuantumComputer.Gate =
        QuantumComputer.gate_extension(gate_add_a, offset, qubit_count - 2)
    subtract_a::QuantumComputer.Gate =
        QuantumComputer.gate_extension(gate_subtract_a, offset, qubit_count - 2)

    cc_add_a::QuantumComputer.Gate =
        QuantumComputer.gate_multi_control(add_a, 2, qubit_count)
    cc_subtract_a::QuantumComputer.Gate =
        QuantumComputer.gate_multi_control(subtract_a, 2, qubit_count)

    cache_path = no_cache ? Array{String,1}(undef, 0) : ["constant_adder_core", string(n_qubit_count + 1), string(n)]
    rebuild && !no_cache && gate_remove_from_cache(cache_path)
    gate_add_n = no_cache ? nothing : QuantumComputer.gate_load_from_cache(cache_path)
    if typeof(gate_add_n) == Nothing
      circuit_add_n::QuantumComputer.Circuit = constant_adder_core(n, n_qubit_count + 1)
      gate_add_n = QuantumComputer.circuit_convert_to_gate(circuit_add_n, cache_path)
    end
    gate_subtract_n::QuantumComputer.Gate = QuantumComputer.gate_invert(gate_add_n)
    add_n::QuantumComputer.Gate =
        QuantumComputer.gate_multi_control(gate_add_n, 1, n_qubit_count + 2)

    c_add_n::QuantumComputer.Gate = QuantumComputer.gate_extension(
        add_n,
        qubit_count - add_n.superposed_qubits_required + 1,
        qubit_count,
    )
    subtract_n::QuantumComputer.Gate = QuantumComputer.gate_extension(
        gate_subtract_n,
        qubit_count - gate_subtract_n.superposed_qubits_required + 1,
        qubit_count,
    )

    gate_qft::QuantumComputer.Gate =
        QuantumComputer.gate_fourier_transform(n_qubit_count + 1)
    gate_inverse_qft::QuantumComputer.Gate = QuantumComputer.gate_invert(gate_qft)

    offset = qubit_count - gate_qft.superposed_qubits_required + 1
    qft::QuantumComputer.Gate =
        QuantumComputer.gate_extension(gate_qft, offset, qubit_count)
    inverse_qft::QuantumComputer.Gate =
        QuantumComputer.gate_extension(gate_inverse_qft, offset, qubit_count)

    gate_cx::QuantumComputer.Gate = QuantumComputer.gate_cx(offset, offset - 1, qubit_count)
    gate_x::QuantumComputer.Gate =
        QuantumComputer.gate_extension(QuantumComputer.gate_x, offset, qubit_count)

    circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()

    QuantumComputer.add_gate_to_circuit!(circuit, cc_add_a)
    QuantumComputer.add_gate_to_circuit!(circuit, subtract_n)
    QuantumComputer.add_gate_to_circuit!(circuit, inverse_qft)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
    QuantumComputer.add_gate_to_circuit!(circuit, qft)
    QuantumComputer.add_gate_to_circuit!(circuit, c_add_n)
    QuantumComputer.add_gate_to_circuit!(circuit, cc_subtract_a)
    QuantumComputer.add_gate_to_circuit!(circuit, inverse_qft)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_x)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_x)
    QuantumComputer.add_gate_to_circuit!(circuit, qft)
    QuantumComputer.add_gate_to_circuit!(circuit, cc_add_a)

    circuit
end

"""
    shor2n3_controlled_modular_multiplier(n, a)

see the modular multiplier circuit in [this paper](https://arxiv.org/pdf/quant-ph/0205095.pdf).

# Arguments:
- `n`: the modulus
- `a`: the constant to multiply by
- `no_cache`: if true, caching will not be employed when constructing gates
- `rebuild`: if true, this will destroy the relevant cached gates before regenerating them
"""
function shor2n3_controlled_modular_multiplier(n::Int64, a::Int64, no_cache = true, rebuild = false)
    n_qubit_count::Int64 = ceil(log2(n))
    qubit_count::Int64 = 2 * n_qubit_count + 3

    circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()

    gate_qft::QuantumComputer.Gate =
        QuantumComputer.gate_fourier_transform(n_qubit_count + 1)
    gate_inverse_qft::QuantumComputer.Gate = QuantumComputer.gate_invert(gate_qft)

    qft::QuantumComputer.Gate =
        QuantumComputer.gate_extension(gate_qft, n_qubit_count + 3, qubit_count)
    inverse_qft::QuantumComputer.Gate =
        QuantumComputer.gate_extension(gate_inverse_qft, n_qubit_count + 3, qubit_count)

    QuantumComputer.add_gate_to_circuit!(circuit, qft)
    for i = 1:n_qubit_count
        control_qubit = n_qubit_count - i + 2
        if control_qubit != 2
            gate_swap::QuantumComputer.Gate =
                QuantumComputer.gate_swap(control_qubit, 2, qubit_count)
            QuantumComputer.add_gate_to_circuit!(circuit, gate_swap)
        end

        cache_path = no_cache ? Array{String,1}(undef, 0) : ["shor2n3", "controlled_controlled_modular_adder", string(qubit_count), string(n), string((a * 2^(i - 1)) % n)]
        rebuild && !no_cache && gate_remove_from_cache(cache_path)
        modular_adder = no_cache ? nothing : QuantumComputer.gate_load_from_cache(cache_path)
        if typeof(modular_adder) == Nothing
            circuit_modular_adder = shor2n3_controlled_controlled_modular_adder(n, (a * 2^(i - 1)) % n, no_cache, rebuild)
            modular_adder = QuantumComputer.circuit_convert_to_gate(circuit_modular_adder, cache_path)
        end

        QuantumComputer.add_gate_to_circuit!(circuit, modular_adder)
        if control_qubit != 2
            QuantumComputer.add_gate_to_circuit!(circuit, gate_swap)
        end
    end
    QuantumComputer.add_gate_to_circuit!(circuit, inverse_qft)

    circuit
end

"""
    shor2n3_controlled_ua(n, a)

the controlled-Ua gate from [this paper](https://arxiv.org/pdf/quant-ph/0205095.pdf).

# Arguments
- `n`: the modulus
- `a`: the constant to multiply `x` by
- `no_cache`: if true, caching will not be employed when constructing gates
- `rebuild`: if true, this will destroy the relevant cached gates before regenerating them
"""
function shor2n3_controlled_ua(n::Int64, a::Int64, no_cache = true, rebuild = false)
    n_qubit_count::Int64 = ceil(log2(n))
    qubit_count::Int64 = 2 * n_qubit_count + 3

    a_inverse::Int64 = invmod(a, n)

    cache_path = no_cache ? Array{String,1}(undef, 0) : ["shor2n3", "controlled_modular_multiplier", string(qubit_count), string(n), string(a)]
    rebuild && !no_cache && gate_remove_from_cache(cache_path)
    modular_multiplier = no_cache ? nothing : QuantumComputer.gate_load_from_cache(cache_path)
    if typeof(modular_multiplier) == Nothing
        circuit_modular_multiplier::QuantumComputer.Circuit =
            QuantumComputer.Circuits.shor2n3_controlled_modular_multiplier(n, a, no_cache, rebuild)
        modular_multiplier = QuantumComputer.circuit_convert_to_gate(circuit_modular_multiplier, cache_path)
    end

    cache_path = no_cache ? Array{String,1}(undef, 0) : ["shor2n3", "controlled_modular_multiplier", string(qubit_count), string(n), string(a_inverse)]
    rebuild && !no_cache && gate_remove_from_cache(cache_path)
    gate_inverse_modular_divider = no_cache ? nothing : QuantumComputer.gate_load_from_cache(cache_path)
    if typeof(gate_inverse_modular_divider) == Nothing
        circuit_inverse_modular_divider::QuantumComputer.Circuit =
            QuantumComputer.Circuits.shor2n3_controlled_modular_multiplier(n, a_inverse, no_cache, rebuild)
        gate_inverse_modular_divider =
            QuantumComputer.circuit_convert_to_gate(circuit_inverse_modular_divider, cache_path)
    end
    modular_divider::QuantumComputer.Gate =
        QuantumComputer.gate_invert(gate_inverse_modular_divider)

    circuit = QuantumComputer.Circuit()

    QuantumComputer.add_gate_to_circuit!(circuit, modular_multiplier)

    # controlled register swap
    swap_matrix::Matrix{Complex{Float64}} = (1.0 + 0.0im) * I(2^(qubit_count - 1))
    for i = 1:n_qubit_count
        gate_swap::QuantumComputer.Gate =
            QuantumComputer.gate_swap(i, n_qubit_count + i + 2, qubit_count - 1)
        # these should all commute
        swap_matrix *= gate_swap.matrix
    end
    gate_swap_registers::QuantumComputer.Gate = QuantumComputer.Gate(swap_matrix)
    c_swap::QuantumComputer.Gate =
        QuantumComputer.gate_multi_control(gate_swap_registers, 1, qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, c_swap)

    QuantumComputer.add_gate_to_circuit!(circuit, modular_divider)

    circuit
end

function shor2n3_execute_hybrid_logic(
    superposition::QuantumComputer.Superposition,
    classical_register::QuantumComputer.ClassicalRegister,
    gates_inverse_qft::Array{QuantumComputer.Gate},
    gate_ua::QuantumComputer.Gate,
    i::Int64,
)
    qubit_count::Int64 = gate_ua.superposed_qubits_required

    # this happens to be the hadamard matrix
    superposition.state = gates_inverse_qft[1].matrix * superposition.state

    # the core gate
    superposition.state = gate_ua.matrix * superposition.state

    # inverse qft
    if i != 1
        for j = 0:(i-2)
            if (classical_register.value & 2^(classical_register.width - j)) != 0
            # if (classical_register.value & 2^j) != 0
                superposition.state = gates_inverse_qft[i-j].matrix * superposition.state
            end
        end
    end
    superposition.state = gates_inverse_qft[1].matrix * superposition.state

    # measure
    measurement = QuantumComputer.Measurement([i], [1], 1)
    # measurement = QuantumComputer.Measurement([classical_register.width - i + 1], [1], 1)
    QuantumComputer.measure_superposition(superposition, classical_register, measurement)

    # conditional pauli x
    if (classical_register.value & 2^(classical_register.width - i + 1)) != 0
    # if (classical_register.value & 2^(i-1)) != 0
        x::QuantumComputer.Gate =
            QuantumComputer.gate_extension(QuantumComputer.gate_x, 1, qubit_count)
        superposition.state = x.matrix * superposition.state
    end
end

"""
    shor2n3_period_finding(n::Int64, a::Int64)

Beauregard's circuit for finding the period of a^x mod n.

# Arguments
- `n`: the modulus
- `a`: the base
- `no_cache`: if true, caching will not be employed when constructing gates
- `rebuild`: if true, this will destroy the relevant cached gates before regenerating them
"""
function shor2n3_period_finding(n::Int64, a::Int64, no_cache = true, rebuild = false)
    n_qubit_count::Int64 = ceil(log2(n))
    qubit_count::Int64 = 2 * n_qubit_count + 3

    gates_ua::Array{Union{QuantumComputer.Gate,Nothing}} =
        Array{Union{QuantumComputer.Gate,Nothing}}(undef, 2 * n_qubit_count)
    gates_inverse_qft::Array{QuantumComputer.Gate} =
        Array{QuantumComputer.Gate}(undef, 2 * n_qubit_count)
    for i = 1:(2*n_qubit_count)
        if i == 1
            gates_inverse_qft[i] =
                QuantumComputer.gate_extension(QuantumComputer.gate_h, 1, qubit_count)
        else
            gates_inverse_qft[i] = QuantumComputer.gate_extension(
                QuantumComputer.gate_p(-pi / 2^(i - 1)),
                1,
                qubit_count,
            )
        end

        local_a = powermod(a, 2^(i - 1), n)
        cache_path = no_cache ? Array{String,1}(undef, 0) : ["shor2n3", "controlled_ua", string(qubit_count), string(n), string(local_a)]
        rebuild && !no_cache && gate_remove_from_cache(cache_path)
        gates_ua[i] = no_cache ? nothing : QuantumComputer.gate_load_from_cache(cache_path)
        if typeof(gates_ua[i]) == Nothing
            circuit_ua::QuantumComputer.Circuit =
                shor2n3_controlled_ua(n, local_a, no_cache, rebuild)
            gates_ua[i] = QuantumComputer.circuit_convert_to_gate(circuit_ua, cache_path)
        end
    end

    circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()

    for i = 1:(2*n_qubit_count)
        component = QuantumComputer.HybridComponent(
            shor2n3_execute_hybrid_logic,
            [gates_inverse_qft, gates_ua[i], i],
        )
        QuantumComputer.add_hybrid_component_to_circuit!(circuit, component)
    end

    circuit
end

"""
    period_finding_for_11x_mod_15()

a 5 qubit circuit that returns a value related to the period of f(x) = 11^x mod 15. this effectively allows implementation of Shor's algorithm for this specific case.
"""
function period_finding_for_11x_mod_15()
    qubit_count = 5

    circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()
    gate_h::QuantumComputer.Gate = QuantumComputer.Gate(
        kron(((1.0 + 0.0im) * I)(4), QuantumComputer.gate_expansion(QuantumComputer.gate_h, 3).matrix),
    )
    QuantumComputer.add_gate_to_circuit!(circuit, gate_h)
    gate_cx::QuantumComputer.Gate = QuantumComputer.gate_cx(3, 4, qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
    gate_cx = QuantumComputer.gate_cx(3, 5, qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
    gate_h = QuantumComputer.gate_extension(QuantumComputer.gate_h, 4, qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_h)
    gate_cp =
        QuantumComputer.gate_control(QuantumComputer.gate_p(pi / 2), [4], 5, qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_cp)
    gate_h = QuantumComputer.gate_extension(QuantumComputer.gate_h, 5, qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_h)

    # uncompute ancilla
    gate_cx = QuantumComputer.gate_cx(3, 2, qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
    gate_cx = QuantumComputer.gate_cx(3, 1, qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
    # done uncomputing

    gate_cp =
        QuantumComputer.gate_control(QuantumComputer.gate_p(pi / 4), [4], 3, qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_cp)
    gate_cp =
        QuantumComputer.gate_control(QuantumComputer.gate_p(pi / 2), [5], 3, qubit_count)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_cp)

    circuit
end

function identification_oracle(qubit_count::Int64, target_input::Int64)
    circuit = QuantumComputer.Circuit()

    matrix_x = QuantumComputer.gate_x.matrix
    matrix_i::Matrix{Complex{Float64}} = ((1.0 + 0.0im)I)(2)

    matrix = matrix_i

    for i in 0:(qubit_count - 2)
        operator = (target_input & (2^i) != 0) ? matrix_i : matrix_x
        matrix = kron(operator, matrix)
    end

    gate_outer::QuantumComputer.Gate = QuantumComputer.Gate(matrix)
    gate_generalized_toffoli::QuantumComputer.Gate = QuantumComputer.gate_multi_control(QuantumComputer.gate_x, qubit_count - 1, qubit_count)

    QuantumComputer.add_gate_to_circuit!(circuit, gate_outer)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_generalized_toffoli)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_outer)

    circuit
end

function grover_diffusion_operator(qubit_count::Int64, gate_h::QuantumComputer.Gate = nothing)
    circuit = QuantumComputer.Circuit()

    if gate_h == nothing
        matrix_h = QuantumComputer.gate_h.matrix
        matrix = matrix_h
        for i in 2:qubit_count
             matrix = kron(matrix, matrix_h)
        end
        gate_h::QuantumComputer.Gate = QuantumComputer.Gate(matrix)
    end

    matrix_x = QuantumComputer.gate_x.matrix
    matrix = matrix_x
    for i in 2:qubit_count
            matrix = kron(matrix, matrix_x)
    end
    gate_x::QuantumComputer.Gate = QuantumComputer.Gate(matrix)

    matrix::Matrix{Complex{Float64}} = ((1.0im) * I)(2)
    identity::Matrix{Complex{Float64}} = ((1.0 + 0.0im) * I)(2^(qubit_count - 2))
    matrix = kron(matrix, identity)
    matrix = kron(matrix, QuantumComputer.gate_h.matrix)
    gate_inner::QuantumComputer.Gate = QuantumComputer.Gate(matrix)

    gate_generalized_toffoli = QuantumComputer.gate_multi_control(QuantumComputer.gate_x, qubit_count - 1, qubit_count)

    QuantumComputer.add_gate_to_circuit!(circuit, gate_h)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_x)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_inner)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_generalized_toffoli)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_inner)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_x)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_h)

    circuit
end

# https://arxiv.org/pdf/quant-ph/0301079.pdf
function grover(oracle::QuantumComputer.Gate, no_cache = true, rebuild = false)
    qubit_count::Int64 = oracle.superposed_qubits_required - 1

    circuit = QuantumComputer.Circuit()

    small_h::QuantumComputer.Gate = QuantumComputer.gate_expansion(QuantumComputer.gate_h, qubit_count)
    gate_h::QuantumComputer.Gate = QuantumComputer.Gate(kron(small_h.matrix, QuantumComputer.gate_h.matrix))
    QuantumComputer.add_gate_to_circuit!(circuit, gate_h)

    cache_path = no_cache ? Array{String,1}(undef, 0) : ["grover", "diffusor", string(qubit_count)]
    rebuild && !no_cache && gate_remove_from_cache(cache_path)
    gate_diffusor = no_cache ? nothing : QuantumComputer.gate_load_from_cache(cache_path)
    if typeof(gate_diffusor) == Nothing
        circuit_diffusor = grover_diffusion_operator(qubit_count, small_h)
        gate_diffusor = QuantumComputer.circuit_convert_to_gate(circuit_diffusor, cache_path)
    end
    gate_extended_diffusor = QuantumComputer.gate_extension(gate_diffusor, 1, qubit_count + 1)

    circuit_grover = QuantumComputer.Circuit()
    QuantumComputer.add_gate_to_circuit!(circuit_grover, oracle)
    QuantumComputer.add_gate_to_circuit!(circuit_grover, gate_extended_diffusor)
    gate_grover = QuantumComputer.circuit_convert_to_gate(circuit_grover)


    iterations = qubit_count < 3 ? 1 : convert(Int64, ceil(pi * ???(2^qubit_count) / 4))
    for _ in 1:iterations
        QuantumComputer.add_gate_to_circuit!(circuit, gate_grover)
    end

    circuit
end

function cx_balanced_oracle(qubit_count::Int64)
    circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()

    for i in 1:(qubit_count - 1)
        gate_cx = QuantumComputer.gate_cx(i, qubit_count, qubit_count)
        QuantumComputer.add_gate_to_circuit!(circuit, gate_cx)
    end

    circuit
end

function constant_oracle(qubit_count::Int64, result::Bool)
    circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()

    gate = result ? QuantumComputer.gate_extension(QuantumComputer.gate_x, qubit_count, qubit_count) : QuantumComputer.Gate(convert(Matrix{Complex{Float64}}, ((1.0 + 0.0im) * I)(2^qubit_count)))
    QuantumComputer.add_gate_to_circuit!(circuit, gate)

    circuit
end

function deutsch_jozsa(oracle::QuantumComputer.Gate)
    qubit_count::Int64 = oracle.superposed_qubits_required - 1

    circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()

    postfix_h = QuantumComputer.gate_expansion(QuantumComputer.gate_h, qubit_count)
    prefix_h = QuantumComputer.Gate(kron(postfix_h.matrix, QuantumComputer.gate_h.matrix))
    postfix_h_extended = QuantumComputer.gate_extension(postfix_h, 1, qubit_count + 1)

    QuantumComputer.add_gate_to_circuit!(circuit, prefix_h)
    QuantumComputer.add_gate_to_circuit!(circuit, oracle)
    QuantumComputer.add_gate_to_circuit!(circuit, postfix_h_extended)

    circuit
end

function bitwise_product_oracle(secret::Int64, qubit_count::Int64)
    circuit = QuantumComputer.Circuit()

    for i in 1:(qubit_count - 1)
        gate = ((secret & 2^(qubit_count - i - 1)) == 0) ? nothing : QuantumComputer.gate_cx(i, qubit_count, qubit_count)
        gate != nothing && QuantumComputer.add_gate_to_circuit!(circuit, gate)
    end

    if length(circuit.components) == 0
        identity = QuantumComputer.Gate(convert(Matrix{Complex{Float64}}, ((1.0 + 0.0im)I)(2^qubit_count)))
        QuantumComputer.add_gate_to_circuit!(circuit, identity)
    end

    circuit
end

function bernstein_vazirani(oracle)
    qubit_count::Int64 = oracle.superposed_qubits_required - 1

    circuit::QuantumComputer.Circuit = QuantumComputer.Circuit()

    gate_h::QuantumComputer.Gate = QuantumComputer.gate_expansion(QuantumComputer.gate_h, qubit_count)
    gate_prefix::QuantumComputer.Gate = QuantumComputer.Gate(kron(gate_h.matrix, QuantumComputer.gate_z.matrix * QuantumComputer.gate_h.matrix))
    gate_postfix::QuantumComputer.Gate = QuantumComputer.gate_extension(gate_h, 1, qubit_count + 1)

    QuantumComputer.add_gate_to_circuit!(circuit, gate_prefix)
    QuantumComputer.add_gate_to_circuit!(circuit, oracle)
    QuantumComputer.add_gate_to_circuit!(circuit, gate_postfix)

    circuit
end

end # module Circuits

end # module QuantumComputer
