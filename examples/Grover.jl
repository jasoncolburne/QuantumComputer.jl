using Pkg
Pkg.activate(".")
using Printf
push!(LOAD_PATH,"../src/")
using QuantumComputer

# configurable runner is below this function
function grover(qubit_count, target_input, sample_size)
    @printf("searching %d inputs for target of %d using Grover's algorithm...\n", 2^qubit_count, target_input)

    println("building circuit...")
    classical_register = QuantumComputer.ClassicalRegister(qubit_count)
    register = QuantumComputer.Register(qubit_count + 1, 1)
    superposition = QuantumComputer.Superposition(register.qubits)

    circuit_oracle = QuantumComputer.Circuits.identification_oracle(qubit_count + 1, target_input)
    gate_oracle = QuantumComputer.circuit_convert_to_gate(circuit_oracle)
    circuit_grover = QuantumComputer.Circuits.grover(gate_oracle, false)
    measurement = QuantumComputer.Measurement(Array(1:qubit_count), Array(1:qubit_count), sample_size)

    circuit = QuantumComputer.Circuit()
    QuantumComputer.add_subcircuit_to_circuit!(circuit, circuit_grover)
    QuantumComputer.add_measurement_to_circuit!(circuit, measurement)

    @printf("sampling (%d times)...\n", sample_size)
    QuantumComputer.apply_circuit_to_superposition!(superposition, circuit, classical_register)

    @printf("classical_register holds value %d\n", classical_register.value)
    classical_register.value == target_input || throw(ErrorException("Grover has failed!"))
    println("success!")
end

# this one is probably fine.
sample_size = 128

qubit_count = length(ARGS) > 0 ? parse(Int64, ARGS[1]) : 10
if length(ARGS) > 1
    QuantumComputer.gate_set_cache_base(ARGS[2])
end
target_input = length(ARGS) > 2 ? parse(Int64, ARGS[3]) : rand(0:(2^qubit_count-1))


grover(qubit_count, target_input, sample_size)
