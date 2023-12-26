using JuMP
using SCS
using Random
using LinearAlgebra
using Statistics
using StatsBase

include("Correlation approximation function.jl")
solver = optimizer_with_attributes(SCS.Optimizer, "max_iters" => 1e+6, "eps_abs" => 1e-8, "eps_rel" => 1e-6)

n = # size of the matrix
    rng1 = MersenneTwister(56308)
M = rand(rng1, -5:5, n, n)
MS = symmetric_matrix(M, n)
MS = convert(Matrix{Float64}, MS)  # Convert it to Matrix{Float64}
f = matrix_to_vector(MS)

#Calculate the Correlation approximation matrix through the following methods: 
#1. Projection method
result, et = @timed correlation_projection_approximation(MS; iteration=1000_000, norm_tolerance=1e-8)
Projex, Projex_iter = result
#2. SDC method
SDC = solveitcorrelationSDC_approx_nodiag(f::Vector{<:Real}, solver)
#3. SQV method
SQV = solveitcorrelationSQV_approx_nodiag(f::Vector{<:Real}, solver)
#4. SQC method
SQC = solveitcorrelationSQC_approx_nodiag(f::Vector{<:Real}, solver)

#Calculate error of our approximation methods
norm(Projex - MS) / norm(MS)
norm(SDC - MS) / norm(MS)
norm(SQVs - MS) / norm(MS)
norm(SQC - MS) / norm(MS)
