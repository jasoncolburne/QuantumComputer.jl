using Pkg
Pkg.activate(".")
using Printf, StatsBase, UnicodePlots
push!(LOAD_PATH,"../src/")
using QuantumComputer

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
    println("building period finding circuit (this can take a few minutes the first time)...")
    shor = QuantumComputer.Circuits.shor2n3_period_finding(n, a, false)

    @printf("sampling circuit (%d samples)...\n", sample_size)
    samples = Array{Int64,1}(undef, sample_size)
    for i in 1:sample_size
        print(".")
        superposition = QuantumComputer.Superposition(qubits)
        QuantumComputer.apply_circuit_to_superposition!(superposition, shor, classical_register)
        samples[i] = classical_register.value
    end
    println()

    labels = Array{String, 1}(undef, 0);
    vals = Array{Integer, 1}(undef, 0);

    counts = countmap(samples)
    for value in sort([key for key in keys(counts)])
        push!(labels, string(string(value), " |", string(value, base=2, pad=2*n_qubit_count), ">"))
        push!(vals, counts[value])
    end

    println(barplot(labels, vals, xlabel = "samples"))
    println()

    filtered = samples
    if length([k for k in keys(counts)]) > 6
        z = median([count for count in values(counts)]) + 1
        z = √(mean([count for count in values(counts)]))
        @printf("eliminating outliers with counts <= %f\n", z)
        filtered = [sample for sample in samples if counts[sample] > z]
        length(filtered) == 0 && throw(ErrorException("found zero non-trivial factors"))
    end

    println("converting samples to phases and determining probable denominators")
    denominators = [denominator(rationalize(value/(2^(2*n_qubit_count)), tol=1/(2^(2*n_qubit_count-1)))) for value in filtered]

    println("removing odd solutions")
    filtered = [denominator for denominator in denominators if denominator % 2 == 0]
    length(filtered) == 0 && throw(ErrorException("found zero non-trivial factors"))

    labels = Array{String, 1}(undef, 0);
    vals = Array{Integer, 1}(undef, 0);

    counts = countmap(filtered)
    for value in sort([key for key in keys(counts)])
        push!(labels, string(value))
        push!(vals, counts[value])
    end
    println(barplot(labels, vals, xlabel = "samples", ylabel = "r"))
    println()

    t(v) = counts[v]
    candidates = sort([v for v in keys(counts)], by = t, rev = true)

    @printf("trying up to %d candidates\n", n_qubit_count)
    firstn(a, n) = a[intersect(eachindex(a), 1:n)]
    for r in firstn(candidates, n_qubit_count)
        @printf("r = %d\n", r)
        raised_a = convert(BigInt, a)^(r >> 1)
        potential_factors = [gcd(raised_a - 1, n), gcd(raised_a + 1, n)]
        factors = [x for x in potential_factors if !(x in [1, n])]

        length(factors) == 0 && continue
        return first(factors)
    end
    throw(ErrorException("found zero non-trivial factors"))
end

# this one is probably fine. play with n and fixed_a
sample_size = 128

# change these (but beware adding another bit to n will increase circuit complexity and memory required)
n = length(ARGS) > 1 ? parse(Int64, ARGS[1]) : 15
fixed_a = length(ARGS) > 1 ? parse(Int64, ARGS[2]) : 7

try
    println()
    x = shor_beauregard(n, fixed_a, sample_size)
    y = convert(Int64, n / x)

    @printf("factored %d into [%d, %d]!\n", n, x, y)
catch
    @printf("zero non-trivial factors of %d found\n", n)
end
