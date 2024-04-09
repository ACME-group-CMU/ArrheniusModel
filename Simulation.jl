function simulate_deposition(fcoeffs, Ea_values, T_range, num_steps, dt)
    # Initialize existing_layers as a 2D array
    new_layer = [1.0 0.0 0.0 0.0 0.0 0.0]
    compositions_all = [] # Initialize compositions for all temperatures
    A = 0.1 # Arrhenius prefactor 
    gprob = 1.0 # Probability of a new layer being deposited
    first_run = true
    for fcoeff in fcoeffs
        println("1")
        compositions_T = [] # Initialize compositions for a single temperature
        for T in T_range
            existing_layers = [1.0 0.0 0.0 0.0 0.0 0.0] # Initial layer(Amorphous)
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
            #println(size(existing_layers))
            push!(compositions_T, existing_layers)
            println("compositions_T: ", size(compositions_T))
        end
        #compositions_T = permutedims(compositions_T)
        if first_run == true
            compositions_all = compositions_T
            first_run = false
        else
            compositions_all = hcat(compositions_all, compositions_T)
        end
        print("compositions_all",size(compositions_all))
    end
    return compositions_all
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