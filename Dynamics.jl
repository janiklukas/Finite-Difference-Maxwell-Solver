"""
The function `plot_fields` calculates and visualizes the electric/magnetic fields from the current and previous potential field states `state_curr` and `state_prev`.
It is intended for the case of charges moving within the xy-plane.
"""
function plot_fields_xy(N_cells::Int64, N_charges::Int64, state_curr::Vector{Float64}, state_prev::Vector{Float64}, Δt::Float64, state_Q::Matrix{Float64})
    #Transform the state vector into grid form
    field_curr = reshape_vector(N_cells, state_curr)
    field_prev = reshape_vector(N_cells, state_prev)
    
    #Calculate the current electromagnetic field state
    field_E = electric_field(N_cells, field_curr, field_prev, Δt)
    field_B = magnetic_field(N_cells, field_curr)
    #field_S = Poynting_field(N_cells, field_E, field_B)
    
    #Visualize the electromagnetic field state
    p1 = plot_slice(N_cells, N_charges, 3, 0.5, field_E, state_Q, "Electric"; col=1) #xy-plane
    p2 = plot_slice(N_cells, N_charges, 1, 0.5, field_B, state_Q, "Magnetic"; col=4) #yz-plane
    #p3 = plot_slice(N_cells, N_charges, 3, 0.5, field_S, state_Q, "Poynting"; col=7) #xy-plane
    
    plot(p1, p2, layout=(1,2), size=(950,500), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
end

"""
The function `moving_charge_xy` calculates and visualizes the time evolution of the electromagnetic field for a charge trajectory given by `Q_array`.
The PDE matrix `A` and the preconditioner `P_inv` must be calculated beforehand using the functions `create_matrix` and `sparse_approximate_inverse`.
The result is saved as `filename`.gif of length `gif_length` with `fps` frames per second.
"""
function moving_charge_xy(N_cells::Int64, A::SparseMatrixCSC{Float64, Int64}, P_inv::SparseMatrixCSC{Float64, Int64}, N_charges::Int64,
    t_max::Float64, steps::Int64, Q_array::Array{Float64, 3}, tol::Float64, maxiters::Int64, gif_length::Float64, fps::Float64, filename::String)

    t_array = collect(range(start=0.0, stop=t_max, length=steps))
    Δt = t_array[2] - t_array[1]
    
    #Determine initial field state
    t0 = time()
    state_prev = initial_state(N_cells, A, P_inv, N_charges, Q_array[1, :, :], Q_array[2, :, :], Δt, tol, maxiters)
    state_curr = copy(state_prev)
    println("Time to compute initial state: ", round(time()-t0, sigdigits=3))
    
    anim = @animate for i in 1:steps-1
        #Plot electric/magnetic fields
        plot_fields_xy(N_cells, N_charges, state_curr, state_prev, Δt, Q_array[i, :, :])
        
        #Calculate the current PDE vector b
        b_curr = create_vector(N_cells, N_charges, Q_array[i, :, :], Q_array[i+1, :, :], Δt)

        #Calculate the next field state
        state_next = Vector(2*state_curr - state_prev + Δt^2 * (A*state_curr + b_curr))
        
        #Shift the previous states
        state_prev = copy(state_curr)
        state_curr = copy(state_next)
    end every round(Int, 1/(gif_length*fps) * length(t_array))

    display(gif(anim, filename * ".gif", fps = fps))
end