function simulate_deposition(compositions_T, fcoeff, Ea_values, T, num_steps, dt)
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
    existing_layers = [1.0 0.0 0.0 0.0 0.0 0.0]
    new_layer = [1.0 0.0 0.0 0.0 0.0 0.0]
    gprob = 1.0 
    #compositions = [] # Initialize compositions for each temperature
    A = 0.1 # Arrhenius prefactor
    k_m = arrhenius_rate_matrix(Ea_values, T) # Calculate rate constants for all phases
    #display(k_m)
    for _ in 1:num_steps
        e = copy(existing_layers)
        for i in 1:size(existing_layers, 1)
            e[i, :] = fcoeff[i] * existing_layers[i, :] # for more inner layers, multiply by smaller fcoeff
        end
        existing_layers += A * e * k_m * dt
        if rand() < gprob # Check if a new layer is deposited
            existing_layers = vcat(reshape(new_layer, 1, :), existing_layers)
            #print(size(existing_layers))
        end
    end
    push!(compositions_T, existing_layers)
    println("compositions_T: ", size(compositions_T))
    return compositions_T    
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