using LinearAlgebra
using Base
using DataFrames
import Plots as plt
using Statistics
using CSV
using StatsBase
using Random
using SCS
using JuMP
using StatsPlots
using Distributions

##note that not all packages above we need it here
include("Correlation approximation function.jl")
solver = optimizer_with_attributes(SCS.Optimizer, "max_iters" => 1e+6, "eps_abs" => 1e-4, "eps_rel" => 1e-4)

df = CSV.read("Credit loan eligiblity data.csv", DataFrame)

# Numeric data types
numerics = ["Int16", "Int32", "Int64", "Float16", "Float32", "Float64"]
num = [Symbol("member_id"), Symbol("loan_amnt"), Symbol("funded_amnt"), Symbol("funded_amnt_inv"), Symbol("int_rate"),
    Symbol("annual_inc"), Symbol("dti"), Symbol("delinq_2yrs"), Symbol("inq_last_6mths"), Symbol("open_acc"),
    Symbol("pub_rec"), Symbol("revol_bal"), Symbol("revol_util"), Symbol("total_acc"), Symbol("total_rec_int"),
    Symbol("total_rec_late_fee"), Symbol("recoveries"), Symbol("collection_recovery_fee"),
    Symbol("collections_12_mths_ex_med"), Symbol("acc_now_delinq"), Symbol("tot_coll_amt"),
    Symbol("tot_cur_bal"), Symbol("total_rev_hi_lim"), Symbol("loan_status")]
numeric_df = select(df, num)
numeric_df = coalesce.(numeric_df, 0)  # Replace missing values with 0
numeric_df = Matrix{Float64}(numeric_df)
correlation_matrix = cor(numeric_df)
#to check the psd of a matrix

# Plot the correlation matrix as a heatmap
s1 = plt.heatmap(correlation_matrix, title="Actual correlation matrix", xticks=(1:length(num), num),
    yticks=(1:length(num), num), xrotation=45, color=:vikO10, clims=(-1, 1))

# add a noise matrix
function add_GUSSIAN_noise(correlation_matrix, noise_level)
    n = size(correlation_matrix, 1)
    #generate a random matrix with gaussian distribution
    rng1 = MersenneTwister(11261)
    random_matrix = randn(rng1, n, n)
    random_matrix = (random_matrix + transpose(random_matrix)) / 2  # Ensure symmetry
    random_matrix[diagind(random_matrix)] .= 0  # Set diagonal elements to zero
    random_matrix_std = std(random_matrix)
    scaled_noise = noise_level * random_matrix / random_matrix_std
    # Scale the noise matrix to have a maximum value of 1 and a minimum value of -1
    scaled_noise = (scaled_noise .- minimum(scaled_noise)) ./ (maximum(scaled_noise) - minimum(scaled_noise)) * 2 .- 1
    # Add the scaled noise to the correlation matrix
    noisy_correlation_matrix = correlation_matrix + scaled_noise
    return noisy_correlation_matrix
end

noise_level = 0.1 # Example noise level

noisy_matrix_G = add_GUSSIAN_noise(correlation_matrix, noise_level)

#plot the heatmap of the noise matrix
s2 = plt.heatmap(noisy_matrix_G, title="Correlation matrix with noise",
    xticks=(1:length(num), num), yticks=(1:length(num), num), xrotation=45, color=:vikO10, clims=(-1, 1))

f = matrix_to_vector(noisy_matrix_G)
result, et = @timed correlation_projection_approximation(noisy_matrix_G; iteration=1000_000, norm_tolerance=1e-8)
Projex, Projex_iter = result
SDC = solveitcorrelationSDH_approx_nodiag(f::Vector{<:Real}, solver)
SQV = solveitcorrelationSQV_approx_nodiag(f::Vector{<:Real}, solver)
SQC = solveitcorrelationSQH_approx_nodiag(f::Vector{<:Real}, solver)

s3 = plt.heatmap(Projex, title="Correlation matrix by Projection method", xticks=(1:length(num), num), yticks=(1:length(num), num), xrotation=45, color=:vikO10, clims=(-1, 1))
s4 = plt.heatmap(SDC, title="Correlation matrix by SDC method", xticks=(1:length(num), num), yticks=(1:length(num), num), xrotation=45, c=:vikO10, clims=(-1, 1))
s5 = plt.heatmap(SQV, title="Correlation matrix by SQV method", xticks=(1:length(num), num), yticks=(1:length(num), num), xrotation=45, c=:vikO10, clims=(-1, 1))
s6 = plt.heatmap(SQC, title="Correlation matrix by SQC method", xticks=(1:length(num), num), yticks=(1:length(num), num), xrotation=45, c=:vikO10, clims=(-1, 1))
plt.plot!(s1, s2, s3, s4, s5, s6, xtickfont=(4), ytickfont=(4), titlefont=("Arial", 8, :bold), layout=(3, 2), size=(1000, 800), title=["" "" "" "" "" ""])

f1 = plt.plot!(title=["Actual correlation matrix " "Correlation matrix with noise" "Correlation matrix by Projection method" "Correlation matrix by SDC method" "Correlation matrix by SQV method" "Correlation matrix by SQC method"])

plt.savefig("vikO1011261.png")


# #error
norm(Projex - correlation_matrix) / norm(correlation_matrix)
norm(SDC - correlation_matrix) / norm(correlation_matrix)
norm(SQV - correlation_matrix) / norm(correlation_matrix)
norm(SQC - correlation_matrix) / norm(correlation_matrix)

norm(Projex - noisy_matrix_G) / norm(noisy_matrix_G)
norm(SDC - noisy_matrix_G) / norm(noisy_matrix_G)
norm(SQV - noisy_matrix_G) / norm(noisy_matrix_G)
norm(SQC - noisy_matrix_G) / norm(noisy_matrix_G)

