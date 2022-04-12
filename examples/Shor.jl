push!(LOAD_PATH,"../src/")

using Printf, QuantumComputer, StatsBase, UnicodePlots

function shor_11x_mod_15(sample_size)
    n_qubit_count = 3
    qubit_count = 5
    n = 15
    a = 11

    @printf("factoring %d using a = %d (optimized %d-qubit circuit)\n", n, a, qubit_count)

    register = QuantumComputer.Register(qubit_count)
    classical_register = QuantumComputer.ClassicalRegister(qubit_count)

    println("building period finding circuit...")
    shor = QuantumComputer.Circuits.period_finding_for_11x_mod_15()
    circuit = QuantumComputer.Circuit()
    measurement = QuantumComputer.Measurement(Array(1:qubit_count), Array(1:qubit_count))

    QuantumComputer.add_subcircuit_to_circuit!(circuit, shor)
    QuantumComputer.add_measurement_to_circuit!(circuit, measurement)

    @printf("sampling circuit (%d samples)...\n", sample_size)
    samples = Array{Int64,1}(undef, sample_size)
    for i in 1:sample_size
        superposition = QuantumComputer.Superposition(register.qubits)
        QuantumComputer.apply_circuit_to_superposition!(superposition, circuit, classical_register)
        samples[i] = classical_register.value & 0x7
    end

    labels = Array{String, 1}(undef, 0);
    values = Array{Integer, 1}(undef, 0);

    for (value, count) in countmap(samples)
      push!(labels, string(string(value), " |", string(value, base=2, pad=n_qubit_count), ">"))
      push!(values, count)
    end

    println(barplot(labels, values, xlabel = "samples"))
    println()

    println("converting samples to phases and determining probable denominators")
    denominators = [denominator(rationalize(value/(2^(n_qubit_count)), tol=1/(2^(2*n_qubit_count-1)))) for value in samples]
    filtered = [denominator for denominator in denominators if !(denominator in [1, n])]

    labels = Array{String, 1}(undef, 0);
    values = Array{Integer, 1}(undef, 0);

    counts = countmap(filtered)
    r = first(counts)[1]
    max_value = first(counts)[1]
    max_count = first(counts)[2]
    for (value, count) in counts
        push!(labels, string(value))
        push!(values, count)

        if count > max_count
            r = value
            max_count = count
        end
        if value > max_value
            max_value = value
        end
    end
    println(barplot(labels, values, xlabel = "samples", ylabel = "r"))
    println()

    r % 2 == 1 && throw(ErrorException("r is odd, quantum algorithm has failed"))
    @printf("found a potential value (%d) for r\n", r)    

    raised_a = a^(r >> 1)
    potential_factors = [gcd(raised_a - 1, n), gcd(raised_a + 1, n)]
    factors = [x for x in potential_factors if !(x in [1, n])]

    length(factors) == 0 && throw(ErrorException("found zero non-trivial factors"))
    first(factors)
end

# configurable runner is below this function
function shor_beauregard(n, a, sample_size)
    n_qubit_count = convert(Int64, ceil(log2(n)))
    @printf("factoring %d using fixed a = %d (generalized %d-qubit circuit)\n", n, fixed_a, 2 * n_qubit_count + 3)
    raw_qubits = [1 0]
    for i in 2:n_qubit_count
      raw_qubits = vcat(raw_qubits, [1 0])
    end
    # here we set x to 1
    raw_qubits = vcat(raw_qubits, [0 1])
    for i in 1:(n_qubit_count + 2)
      raw_qubits = vcat(raw_qubits, [1 0])
    end

    qubits = convert(Array{Complex{Float64},2}, raw_qubits)
    classical_register = QuantumComputer.ClassicalRegister(2 * n_qubit_count)
    println("building period finding circuit (this can take a few minutes)...")
    shor = QuantumComputer.Circuits.shor2n3_period_finding(n, a)

    @printf("sampling circuit (%d samples)...\n", sample_size)
    samples = Array{Int64,1}(undef, sample_size)
    for i in 1:sample_size
        superposition = QuantumComputer.Superposition(qubits)
        QuantumComputer.apply_circuit_to_superposition!(superposition, shor, classical_register)
        samples[i] = classical_register.value
    end

    labels = Array{String, 1}(undef, 0);
    values = Array{Integer, 1}(undef, 0);

    for (value, count) in countmap(samples)
      push!(labels, string(string(value), " |", string(value, base=2, pad=2*n_qubit_count), ">"))
      push!(values, count)
    end

    println(barplot(labels, values, xlabel = "samples"))
    println()

    println("converting samples to phases and determining probable denominators")
    denominators = [denominator(rationalize(value/(2^(2*n_qubit_count)), tol=1/(2^(2*n_qubit_count-1)))) for value in samples]

    println("removing trivial solutions")
    filtered = [denominator for denominator in denominators if !(denominator in [1, n])]

    println("eliminating outliers")
    counts = countmap(filtered)
    max_count = first(counts)[2]
    for (value, count) in counts
        if count > max_count
            max_count = count
        end
    end
    filtered = [sample for sample in filtered if counts[sample] >= 2*n_qubit_count - 1]

    labels = Array{String, 1}(undef, 0);
    values = Array{Integer, 1}(undef, 0);

    counts = countmap(filtered)
    r = first(counts)[1]
    max_value = first(counts)[1]
    max_count = first(counts)[2]
    for (value, count) in counts
        push!(labels, string(value))
        push!(values, count)

        if count > max_count
            r = value
            max_count = count
        end
        if value > max_value
            max_value = value
        end
    end
    println(barplot(labels, values, xlabel = "samples", ylabel = "r"))
    println()

    r % 2 == 1 && throw(ErrorException("r is odd, quantum algorithm has failed"))
    @printf("found a potential value (%d) for r\n", r)

    raised_a = a^(r >> 1)
    potential_factors = [gcd(raised_a - 1, n), gcd(raised_a + 1, n)]
    factors = [x for x in potential_factors if !(x in [1, n])]

    length(factors) == 0 && throw(ErrorException("found zero non-trivial factors"))
    first(factors)
end


# this one is probably fine. play with n and fixed_a
sample_size = 128

x = shor_11x_mod_15(sample_size)
y = convert(Int64, 15 / x)

@printf("factored 15 into [%d, %d]!\n", x, y)

println()

# change these (but beware adding another bit to n will increase circuit complexity and memory required)
n = 15
fixed_a = 7

x = shor_beauregard(n, fixed_a, sample_size)
y = convert(Int64, n / x)

@printf("factored %d into [%d, %d]!\n", n, x, y)