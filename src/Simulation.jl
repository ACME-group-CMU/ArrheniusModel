function meshgrid(x, y)
"""
To build the meshgrids for the x and y values (Temperature and Flow Rate)
"""
    x_grid = repeat(reshape(x, 1, :), length(y), 1)
    y_grid = repeat(y, 1, length(x))
    return x_grid, y_grid
end

function deposition_rates!(dc, c, p, t)
    # Unpack parameters
    fcoeff = p.fcoeff
    K = p.K
    j0 = p.j0
    j = p.j
    dt = p.dt
    num_steps = Int(p.num_steps)
    num_layers = Int(p.num_layers)
    # Calculate deposition rates
    j = floor(Int, t / 0.5) + 1
    f = reverse(fcoeff[j: num_layers+j-1])
    dc .= c .* f * K
    if j != j0
        c[j+1, 1] = 1.0
        j = j0
    end
end

function simulate_deposition(flow_rate, T, barriers::Matrix, para_sim, decay_constant = 0.00001, final_step=true)
    num_steps, num_layers, dt = para_sim
    # Initialize existing_layers as a 2D array
    decay_coefficients = decay_constant * flow_rate
    fcoeff = flow_coefficient("exponential", num_layers, decay_coefficients)
    n = size(barriers, 1)
    c0 = zeros(num_layers, n)
    c0[1, 1] = 1.0
    K = arrhenius_rate(barriers, T)
    j = 0
    j0 = 0
    p = (fcoeff, K, j0, j, dt, num_steps, num_layers)
    p = NamedTuple{(:fcoeff, :K, :j0, :j, :dt, :num_steps, :num_layers)}(p)
    p = ComponentArray(p)
    tspan = (0.0, (num_steps-1) * dt)
    prob = ODEProblem(deposition_rates!, c0, tspan, p)
    if final_step
        sol = solve(prob, Euler(), save_everystep = false, dt = dt)
        return sol.u[end]
    else
        sol = solve(prob, Euler(), saveat = 0.5, dt = dt)
        return sol.u
    end
end

function simulate_deposition!(sol::Matrix, flow_rate::Vector, T::Vector, barriers::Matrix, para_sim)
    sol .= simulate_deposition(flow_rate, T, barriers, para_sim)
    return nothing
end

function oldsimulate_deposition(fcoeff, pe::PhaseEnergies, T, num_steps, dt)
"""
Inputs:
    fcoeff: flow coefficients for each phase
    pe: Phase Energires datastructure and saved values
Output:
    Layers of Phases according to one temperature and flow-rate conditions
"""
    # Initialize existing_layers as a 2D array
    n = n_phases(pe)
    new_layer = zeros(1, n)
    new_layer[1] = 1.0
    existing_layers = new_layer
    gprob = 1.0 
    #compositions = [] # Initialize compositions for each temperature
    A = 0.1 # Arrhenius prefactor
    arrhenius_rate(pe, T) # Calculate rate constants for all phases
    #display(k_m)
    for n in 1:num_steps
        e = copy(existing_layers)
        for i in 1:size(existing_layers, 1)
            e[i, :] = fcoeff[i] * existing_layers[i, :] # for more inner layers, multiply by smaller fcoeff
        end
        existing_layers += A * e * pe.K * dt
        #println("existing_layers: ", existing_layers)  # Print the values of existing_layers
        if rand() < gprob # Check if a new layer is deposited
            existing_layers = vcat(reshape(new_layer, 1, :), existing_layers)
        end
        # Store the composition of the top layer at this time step
    end
    return existing_layers
end



function most_preferable_state(composition, threshold, phase_names)
    # Average all the layers
    avg_composition = mean(composition, dims=1)
    # True if greater than threshold
    boolean_composition = avg_composition .> threshold
    # Create a string with the names of the phases where the composition is true
    phase_string = join([phase_names[i] for i in 1:length(boolean_composition) if boolean_composition[i]], "+")
    return phase_string
end