# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3
#     name: python3
# ---

# + id="OSOaG-m77vLy"
function arrhenius_rate(A, Ea, T)
    k = A * exp(-Ea / (8.314 * T))  # R = 8.314 J/(mol*K)
    return k
end

function arrhenius_transform_with_deposition(existing_layers, new_layer, A_values, Ea_values, T)
    new_layers = []

    # Transform based on Arrhenius equation
    for existing_layer in existing_layers
        # Calculate transformation rates for existing layer
        k_values = [arrhenius_rate(A, Ea, T) for (A, Ea) in zip(A_values, Ea_values)]

        # Transform existing layer
        transformed_layer = similar(existing_layer)
        for i in 1:length(existing_layer)
            transformed_layer[i] = existing_layer[i] * (1 - k_values[i])
            for j in 1:length(existing_layer)
                if i != j
                    transformed_layer[i] += existing_layer[j] * k_values[j]
                end
            end
        end
        push!(new_layers, transformed_layer)
    end

    # Deposit new layer
    push!(new_layers, new_layer)

    return new_layers
end
