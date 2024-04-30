function meshgrid(x, y)
    x_grid = repeat(reshape(x, 1, :), length(y), 1)
    y_grid = repeat(y, 1, length(x))
    return x_grid, y_grid
end

function simulate_deposition(fcoeff, pe::PhaseEnergies, T, num_steps, dt)
"""
Inputs:
    Compositions_T: empty 3D array of compositions for each temperature
                    (n layers of composition of each phase [x, a ,b ,k, y, d])
    fcoeff: flow coefficients for each phase
    Ea_values: Activation energies for each phase
Output:
    Compositions_T: 3D array of compositions for each temperature
"""
    # Initialize existing_layers as a 2D array
    n = n_phases(pe)
    existing_layers = [1.0 0.0 0.0 0.0 0.0 0.0]
    new_layer = [1.0 0.0 0.0 0.0 0.0 0.0]
    gprob = 1.0 
    #compositions = [] # Initialize compositions for each temperature
    A = 0.1 # Arrhenius prefactor
    k_m = arrhenius_rate(pe, T) # Calculate rate constants for all phases
    #display(k_m)
    for n in 1:num_steps
        e = copy(existing_layers)
        for i in 1:size(existing_layers, 1)
            e[i, :] = fcoeff[i] * existing_layers[i, :] # for more inner layers, multiply by smaller fcoeff
        end
        existing_layers += A * e * k_m * dt
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