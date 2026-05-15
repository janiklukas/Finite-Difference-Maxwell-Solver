"""
The function `electric_field` computes the electric field from the current and previous field states `field_curr` and `field_prev`.
"""
function electric_field(N_cells::Int64, field_curr::Array{Float64, 4}, field_prev::Array{Float64, 4}, Δt::Float64)
    E = zeros(N_cells, N_cells, N_cells, 3)
    ϕ = field_curr[:, :, :, 1]
    a = 1.0/N_cells
    
    #∇ϕ term
    for d in 1:3
        E[:, :, :, d] = - (selectdim(ϕ, d, vcat(N_cells:N_cells, 1:N_cells-1)) - selectdim(ϕ, d, vcat(2:N_cells, 1:1)))/a
    end

    #∂tA term
    E[:, :, :, 1:3] += - (field_curr[:, :, :, 2:4] - field_prev[:, :, :, 2:4])/Δt

    return E
end

"""
The function `magnetic_field` computes the electric field from the current field state `field_curr`.
"""
function magnetic_field(N_cells::Int64, field_curr::Array{Float64, 4})
    B = zeros(N_cells, N_cells, N_cells, 3)
    A_x = field_curr[:, :, :, 2]
    A_y = field_curr[:, :, :, 3]
    A_z = field_curr[:, :, :, 4]
    a = 1.0/N_cells
    
    #x component
    B[:, :, :, 1] += (selectdim(A_z, 2, vcat(N_cells:N_cells, 1:N_cells-1)) - selectdim(A_z, 2, vcat(2:N_cells, 1:1)))/a
    B[:, :, :, 1] += - (selectdim(A_y, 3, vcat(N_cells:N_cells, 1:N_cells-1)) - selectdim(A_y, 3, vcat(2:N_cells, 1:1)))/a

    #y component
    B[:, :, :, 2] += (selectdim(A_x, 3, vcat(N_cells:N_cells, 1:N_cells-1)) - selectdim(A_x, 3, vcat(2:N_cells, 1:1)))/a
    B[:, :, :, 2] += - (selectdim(A_z, 1, vcat(N_cells:N_cells, 1:N_cells-1)) - selectdim(A_z, 1, vcat(2:N_cells, 1:1)))/a
    
    #z component
    B[:, :, :, 3] += (selectdim(A_y, 1, vcat(N_cells:N_cells, 1:N_cells-1)) - selectdim(A_y, 1, vcat(2:N_cells, 1:1)))/a
    B[:, :, :, 3] += - (selectdim(A_x, 2, vcat(N_cells:N_cells, 1:N_cells-1)) - selectdim(A_x, 2, vcat(2:N_cells, 1:1)))/a

    return B
end

"""
The function `Poynting_field` computes the Poynting field from the current electric/magnetic field state `field_E` and `field_B`.
"""
function Poynting_field(N_cells::Int64, field_E::Array{Float64, 4}, field_B::Array{Float64, 4})
    S = zeros(N_cells, N_cells, N_cells, 3)

    #x component
    S[:, :, :, 1] = @. field_E[:, :, :, 2] * field_B[:, :, :, 3] - field_E[:, :, :, 3] * field_B[:, :, :, 2]

    #y component
    S[:, :, :, 2] = @. field_E[:, :, :, 3] * field_B[:, :, :, 1] - field_E[:, :, :, 1] * field_B[:, :, :, 3]

    #z component
    S[:, :, :, 3] = @. field_E[:, :, :, 1] * field_B[:, :, :, 2] - field_E[:, :, :, 2] * field_B[:, :, :, 1]

    return S
end

"""
The function `plot_slice` visualizes a 2D slice of a given `field` as well as the charges positions in `state_Q`.
The parameter `fixed` determines which coordinate is kept constant at the value `d_0`.
"""
function plot_slice(N_cells::Int64, N_charges::Int64, fixed::Int64, d_0::Float64, field::Array{Float64, 4}, state_Q::Matrix{Float64}, name::String; col::Int64=1)
    @assert fixed in [1, 2, 3] ""
    a = 1.0/N_cells
    k = round(Int, d_0 / a)
    
    #Calculate 2D field slice
    if fixed == 1
        str = " (x/L = "
        slice = field[k, :, :, 1:2]
    elseif fixed == 2
        str = " (y/L = "
        slice = field[:, k, :, 1:2]
    elseif fixed == 3
        str = " (z/L = "
        slice = field[:, :, k, 1:2]
    end
    
    #Normalize the field vectors
    norms = [norm(slice[i, j, :]) for i in 1:N_cells for j in 1:N_cells]
    slice = 0.8 * a * slice / maximum(norms)
    
    #Initialize plot
    p = plot(legend=false, title=name*str*string(round(a*k, sigdigits = 2))*")", aspect_ratio=:equal)
    plot!(xticks=[0, 1.0], yticks=[0, 1.0], minorticks=N_cells, minorgrid=true)
    plot!(size=(400,400), left_margin=-2Plots.mm, bottom_margin=-2Plots.mm)
    xlims!(0, 1.0); ylims!(0, 1.0)
    if fixed == 1
        xlabel!(L"y/L"); ylabel!(L"z/L")
    elseif fixed == 2
        xlabel!(L"x/L"); ylabel!(L"z/L")
    elseif fixed == 3
        xlabel!(L"x/L"); ylabel!(L"y/L")
    end
    
    #Plot field for each cell
    for i in 1:N_cells
        for j in 1:N_cells
            P = a * [i-1/2, j-1/2] #cell center location
            Q = P + slice[i, j, :]
            plot!([P[1], Q[1]], [P[2], Q[2]], color=col, linewidth=2, arrow=:closed)
        end
    end

    #Add charges (semi-transparent if outside plot plane)
    for n in 1:N_charges
        if round(Int, state_Q[n, 1+fixed] / a) == k
            alpha = 1.0
        else
            alpha = 0.5
        end
        coords = filter(x -> x != fixed, [1, 2, 3])
        if state_Q[n, 1] > 0
            scatter!([state_Q[n, 1+coords[1]]], [state_Q[n, 1+coords[2]]], color=2, alpha=alpha)
        else
            scatter!([state_Q[n, 1+coords[1]]], [state_Q[n, 1+coords[2]]], color=1, alpha=alpha)
        end
    end

    return p
end