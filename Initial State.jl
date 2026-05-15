"""
The function `create_matrix` calculates the sparse symmetric PDE matrix `A` of size 4*`N_cells`^3.
"""
function create_matrix(N_cells::Int64)
    t0 = time()

    #Initialize empty matrix
    n = 4*N_cells^3
    matrix = spzeros(n, n)
    
    #Save directions for neighbor calculations
    directions = [CartesianIndex(1, 0, 0), CartesianIndex(0, 1, 0), CartesianIndex(0, 0, 1)]
    signs = [-1, +1]
    
    #Compute the i-th matrix column
    for i in 1:n
        e_i = zeros(n)
        e_i[i] = 1
        
        #Find the nonzero entry in the corresponding field state vector
        field_state = reshape(e_i, (N_cells, N_cells, N_cells, 4)) #not needed explicitly!
        nonzero = findall(x -> x != 0, field_state)[1]
        
        #Save relevant info for neighbor calculations
        cell = CartesianIndex(nonzero[1], nonzero[2], nonzero[3]) #cell location
        comp = nonzero[4] #field component
        eqs = zeros(N_cells, N_cells, N_cells, 4) #4 eqs per cell

        #Go through neighboring cells and save contributions to eqs
        for d in 1:3
            for σ in signs
                neighbor = cell + σ * directions[d]
                #inner cells:
                if neighbor[d] > 0 && neighbor[d] < N_cells + 1
                    eqs[neighbor, comp] = 1
                #boundary cells (periodic BCs):
                elseif neighbor[d] == 0
                    neighbor = neighbor + N_cells * directions[d]
                    eqs[neighbor, comp] = 1
                elseif neighbor[d] == N_cells + 1
                    neighbor = neighbor - N_cells * directions[d]
                    eqs[neighbor, comp] = 1
                end
            end
        end
        
        #Add the final Laplacian term
        eqs[cell, comp] = -6
        
        #Overwrite the corresponding matrix column
        col = sparse(reshape(eqs, n))
        inds = findnz(col)[1]
        for ind in inds
            matrix[ind, i] = col[ind]
        end
    end

    println("Time to compute A: ", round(time()-t0, sigdigits=3))

    return matrix
end

"""
The function `create_vector` calculates the sparse PDE vector `b` of length 4*`N_cells`^3 from the current and next charge state `Q_curr` and `Q_next`.
"""
function create_vector(N_cells::Int64, N_charges::Int64, Q_curr::Matrix{Float64}, Q_next::Matrix{Float64}, Δt::Float64)
    #Initialize equation vector
    n = 4*N_cells^3
    eqs = zeros(N_cells, N_cells, N_cells, 4)
    
    #Determine charge locations on the grid
    a = 1.0/N_cells
    cells = ceil.(Int, Q_curr[:, 2:4] / a)

    #Compute total charge & current for every grid cell
    for n in 1:N_charges
        eqs[cells[n, :]..., 1] += 4π / a^3 * Q_curr[n, 1]
        eqs[cells[n, :]..., 2:4] += 4π / a^3 * Q_curr[n, 1] * (Q_next[2:4] - Q_curr[2:4])/Δt
    end
    
    return sparse(reshape(eqs, n))
end

"""
The function `sparse_ansatz` produces a sparse symmetric matrix of size 4*`N_cells`^3 whose diagonals are filled with the values in `params`.
"""
function sparse_ansatz(N_cells::Int64, params::Vector{Float64})
    n = 4*N_cells^3
    p = length(params)
    @assert p >= 3 "The number of `params` must be at least 3!"
    
    #Populate central diagonals
    pairs = [0 => fill(params[1], n)]
    for k in 1:p-2
        push!(pairs, k => fill(params[k+1], n-k))
        push!(pairs, -k => fill(params[k+1], n-k))
    end
    
    #Populate boundary diagonal
    push!(pairs, N_cells => fill(params[end], n-N_cells))
    push!(pairs, -N_cells => fill(params[end], n-N_cells))

    P_inv = spdiagm(n, n, pairs...)

    return P_inv
end

"""
The function `sparse_approximate_inverse` computes an approximate inverse `P_inv` of the PDE matrix `A` based on the function `sparse_ansatz`.
"""
function sparse_approximate_inverse(A::SparseMatrixCSC{Float64, Int64}, N_cells::Int64, nparams::Int64)
    t0 = time()
    n = 4*N_cells^3

    #Determine optimal sparse inverse
    function f(params)
        M = sparse_ansatz(N_cells, params)*A - sparse(1:n, 1:n, ones(n))

        return sqrt(abs(tr(M'*M)))
    end
    res = optimize(f, zeros(nparams))
    params = Optim.minimizer(res)
    P_inv = sparse_ansatz(N_cells, params)

    println("Time to compute P_inv: ", round(time()-t0, sigdigits=3))

    return P_inv
end

"""
The function `conjugate_residual` solves the linear system Ax=b for x, where `A` is a symmetric sparse matrix.
To ensure stability, a preconditioner `P_inv` should be supplied. The algorithm terminates when |Ax-b| <= tol*|b| or after `maxiters` iterations.
"""
function conjugate_residual(A::SparseMatrixCSC{Float64, Int64}, b::SparseVector{Float64, Int64}, P_inv::SparseMatrixCSC{Float64, Int64}, tol::Float64, maxiters::Int64)
    @assert norm(A-transpose(A)) < 1e-1 "The matrix `A` must be symmetric!"
    
    #Initialize the relevant variables
    x = spzeros(length(b))
    r = P_inv * (b - A*x)
    p = copy(r)
    Ax = A*x
    Ar = A*r
    Ap = copy(Ar)
    iters = 0

    #Iterate until the desired accuracy is achieved
    for _ in 1:maxiters
        iters += 1
        α = dot(r, Ar) / dot(Ap, P_inv * Ap)
        β = 1 / dot(r, Ar)
        x += α * p
        Ax += α * Ap
        if norm(Ax - b) <= tol * norm(b)
            break
        end
        r -= α * P_inv * Ap
        Ar = A*r
        β *= dot(r, Ar)
        p = r + β * p
        Ap = Ar + β * Ap
    end

    if iters == maxiters error("Maximum iterations reached!") end

    return x
end

"""
The function `reshape_vector` converts a vector of length 4*`N_cells`^3 to a 4-dimensional array that corresponds to the spatial arrangement of the grid cells.
"""
function reshape_vector(N_cells::Int64, vec::Vector{Float64})
    return reshape(vec, (N_cells, N_cells, N_cells, 4))
end

"""
The function `initial_state` solves the static Maxwell potential equations for a charge state given by `Q_init` and `Q_next`. 
The PDE matrix `A` and the preconditioner `P_inv` must be calculated beforehand using the functions `create_matrix` and `sparse_approximate_inverse`.
"""
function initial_state(N_cells::Int64, A::SparseMatrixCSC{Float64, Int64}, P_inv::SparseMatrixCSC{Float64, Int64}, N_charges::Int64, Q_init::Matrix{Float64}, Q_next::Matrix{Float64}, Δt::Float64, tol::Float64, maxiters::Int64)
    #Calculate the PDE vector b
    b = create_vector(N_cells, N_charges, Q_init, Q_next, Δt)
    
    #Solve the static field equations
    state_init = conjugate_residual(-A, b, P_inv, tol, maxiters)
    
    return Vector(state_init)
end