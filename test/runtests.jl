import QuantumComputer
using Test

@testset "Register" begin
  register = QuantumComputer.Register(1)
  @test register.qubits == [1 0]

  register = QuantumComputer.Register(1, 1)
  @test register.qubits == [0 1]

  register = QuantumComputer.Register(2, 2)
  @test register.qubits == [0 1; 1 0]

  register = QuantumComputer.Register(10)
  @test size(register.qubits) == (10, 2)
end

@testset "Superposition" begin
  register = QuantumComputer.Register(3, 7)
  superposition = QuantumComputer.Superposition(register.qubits)
  @test length(superposition.state) == 2^3
  @test superposition.state[8] == 1
  @test sum(superposition.state) == 1
end

@testset "Measurement" begin
  @test_throws DomainError QuantumComputer.Measurement([1, 1], [1])
  @test_throws DomainError QuantumComputer.Measurement([1], [1, 1])
  @test_throws DomainError QuantumComputer.Measurement([1, 1], [1, 2])
  @test_throws DomainError QuantumComputer.Measurement([1, 2], [1, 1])
end

@testset "qubits_operated_on_by_unitary_matrix" begin
  identity_matrix = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
  @test QuantumComputer.qubits_operated_on_by_unitary_matrix(identity_matrix) == 2
  @test_throws DomainError QuantumComputer.qubits_operated_on_by_unitary_matrix(identity_matrix[1:3, :])
  @test_throws DomainError QuantumComputer.qubits_operated_on_by_unitary_matrix(identity_matrix[:, 1:3])
end

@testset "Gate" begin
  identity_matrix = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
  gate::QuantumComputer.Gate = QuantumComputer.Gate(identity_matrix)
  @test gate.matrix == identity_matrix
  @test gate.superposed_qubits_required == 2
end

@testset "gate_x" begin
  state::Array{Complex{Float64}, 1} = [1, 0]
  @test QuantumComputer.gate_x.matrix * state == [0, 1]

  state = [0, 1]
  @test QuantumComputer.gate_x.matrix * state == [1, 0]
end

@testset "gate_cx" begin
  gate::QuantumComputer.Gate = QuantumComputer.gate_cx(1, 3, 3)

  qubits::Matrix{Complex{Float64}} = [1 0; 1 0; 1 0]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == state

  qubits = [1 0; 1 0; 0 1]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == state

  qubits = [0 1; 1 0; 1 0]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == [0, 0, 0, 0, 0, 1, 0, 0]

  qubits = [0 1; 1 0; 0 1]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == [0, 0, 0, 0, 1, 0, 0, 0]

  gate = QuantumComputer.gate_cx(1, 2, 3)

  qubits = [0 1; 1 0; 0 1]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == [0, 0, 0, 0, 0, 0, 0, 1]
end

@testset "gate_cnx" begin
  gate::QuantumComputer.Gate = QuantumComputer.gate_cnx([1, 2], 3, 3)

  qubits::Matrix{Complex{Float64}} = [0 1; 1 0; 1 0]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == state

  qubits = [1 0; 0 1; 1 0]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == state

  qubits = [1 0; 1 0; 1 0]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == state

  qubits = [0 1; 0 1; 1 0]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == [0, 0, 0, 0, 0, 0, 0, 1]
end

@testset "gate_swap" begin
  gate = QuantumComputer.gate_swap(1, 3, 3)

  qubits::Matrix{Complex{Float64}} = [1 0; 1 0; 1 0]
  state::Array{Complex{Float64}} = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == state

  qubits = [0 1; 1 0; 0 1]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == state

  qubits = [0 1; 1 0; 1 0]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == [0, 1, 0, 0, 0, 0, 0, 0]

  qubits = [1 0; 1 0; 0 1]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == [0, 0, 0, 0, 1, 0, 0, 0]
end

@testset "gate_multi_control" begin
  inner_gate = QuantumComputer.gate_swap(1, 3, 3)
  gate = QuantumComputer.gate_multi_control(inner_gate, 1, 4)

  qubits::Matrix{Complex{Float64}} = [1 0; 0 1; 1 0; 1 0]
  state::Array{Complex{Float64}} = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == state

  qubits = [1 0; 1 0; 1 0; 0 1]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == state

  qubits = [0 1; 0 1; 1 0; 1 0]
  state = QuantumComputer.qubit_tensor_product(qubits)
  expected_qubits::Matrix{Complex{Float64}} = [0 1; 1 0; 1 0; 0 1]
  expected_state::Array{Complex{Float64}} = QuantumComputer.qubit_tensor_product(expected_qubits)
  @test gate.matrix * state == expected_state

  qubits = [0 1; 1 0; 0 1; 0 1]
  state = QuantumComputer.qubit_tensor_product(qubits)
  expected_qubits = [0 1; 0 1; 0 1; 1 0]
  expected_state = QuantumComputer.qubit_tensor_product(expected_qubits)
  @test gate.matrix * state == expected_state
end

@testset "gate_fourier_transform" begin
  gate = QuantumComputer.gate_fourier_transform(2)

  qubits::Matrix{Complex{Float64}} = [1 0; 1 0]
  state::Array{Complex{Float64}} = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state == [0.5, 0.5, 0.5, 0.5]

  qubits = [0 1; 0 1]
  state = QuantumComputer.qubit_tensor_product(qubits)
  @test gate.matrix * state ≈ [0.5, -0.5im, -0.5, 0.5im] atol=0.000000000000001
end