
function symmetric_matrix(X, n)
    for j in 1:n
        for i in 1:n
            X[i, j] = X[j, i]
        end
    end
    return value.(X)
end

function matrix_to_vector(A::Array{T,2} where {T<:Real})

    n = size(A, 1)
    v = Vector{Float64}(undef, n * (n - 1) ÷ 2)  # Preallocate vector

    k = 1
    for j in 1:n
        for i in 1:j-1
            v[k] = A[i, j]
            k += 1
        end
    end

    return v
end

function project_to_PSD(F)
    n = size(F, 1)
    λ, U = eigen(F)
    λ = real(λ)
    eigenvecs = real.(U)
    λ[λ.<0] .= 0

    S = U * Diagonal(λ) * U'
    return S

end

function correlation_projection_approximation(A; iteration=1000_000, norm_tolerance=1e-8)
    n = size(A, 1)
    ONES = ones(n)
    P_K(X, k) = X - Diagonal(diag(X)) + Diagonal(ONES)
    P_S(X) = project_to_PSD(X)
    A_k = copy(A)
    for k in 1:iteration
        A_prev = A_k
        A_k = A_k + P_K(P_S(A_k), k) - P_S(A_k)
        A_k = A_k - Diagonal(diag(A_k)) + Diagonal(ONES)
        λ, U = eigen(A_k)
        λ = real(λ)
        A_k = whatever(U, λ, ONES)
        # Check for convergence based on norm difference
        norm_diff = norm(A_k - A_prev)
        if norm_diff ≤ norm_tolerance
            return A_k, k
        end
    end

    return A_k, iteration
end

function whatever(U, λ, ONES)
    U = real(U)
    λ[λ.<0] .= 0
    new_A = U * Diagonal(λ) * U'
    new_A = new_A - Diagonal(diag(new_A)) + Diagonal(ONES)
    return new_A
end

function solveitcorrelationSDC_approx_nodiag(f::Vector{<:Real}, solver)
    # Finding the approximation correlation matrix by using the m by 1 vector f through SDH method, where m=n(n-1)/2
    m, = size(f)
    n = Int(round((0.5 + sqrt(1 + 8 * m)) / 2))
    model = Model(solver)

    @variable(model, t >= 0)
    @variable(model, v[1:m])
    @variable(model, X[1:n, 1:n], PSD)

    @constraint(model, [i in 1:n], X[i, i] == 1)
    triu_elements = [X[i, j] for j in 2:n for i in 1:(j-1)]
    @constraint(model, svecC, v == [sqrt(2) * val for (i, val) in enumerate(f)] .* (f .- triu_elements))
    @constraint(model, [I(m) v; v' t] in PSDCone())

    @objective(model, Min, t)
    optimize!(model)

    return value.(X)
    println(solution_summary(model))
end

function solveitcorrelationSQC_approx_nodiag(f::Vector{<:Real}, solver)
    # Finding the approximation correlation matrix by using the m by 1 vector f through SQH method, where m=n(n-1)/2
    m, = size(f)
    n = Int(round((0.5 + sqrt(1 + 8 * m)) / 2))
    model = Model(solver)

    @variable(model, t >= 0)
    @variable(model, v[1:m])
    @variable(model, X[1:n, 1:n], PSD)

    @constraint(model, [i in 1:n], X[i, i] == 1)
    triu_elements = [X[i, j] for j in 2:n for i in 1:(j-1)]
    @constraint(model, svecC, v == [sqrt(2) * val for (i, val) in enumerate(f)] .* (f .- triu_elements))
    @constraint(model, [t; v] in SecondOrderCone())

    @objective(model, Min, t)
    optimize!(model)

    return value.(X)
    println(solution_summary(model))
end

function solveitcorrelationSQV_approx_nodiag(f::Vector{<:Real}, solver)
    # Finding the approximation correlation matrix by using the m by 1 vector f through SQV method, where m=n(n-1)/2
    m, = size(f)
    n = Int(round((0.5 + sqrt(1 + 8 * m)) / 2))
    model = Model(solver)

    @variable(model, t >= 0)
    @variable(model, X[1:n, 1:n], PSD)

    @constraint(model, [i in 1:n], X[i, i] == 1)
    triu_elements = [X[i, j] for j in 2:n for i in 1:(j-1)]
    @constraint(model, [t; f .- triu_elements] in SecondOrderCone())

    @objective(model, Min, t)
    println(model)
    optimize!(model)

    return value.(X)
    println(solution_summary(model))
end

