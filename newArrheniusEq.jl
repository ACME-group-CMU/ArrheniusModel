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
#include("Energies.jl")

kb = 8.617e-5 #eV/K
function arrhenius_rate(A, Ea, T)
    k = A * exp(-Ea / (kb* T))
    return k
end

function arrhenius_rate_matrix(A, Ea_m, T, dt)
    k_matrix = [arrhenius_rate(A, Ea_m[i, j], T) for i in 1:size(Ea_m, 1), j in 1:size(Ea_m, 2)]
    k = copy(k_matrix)
    # Adjust the diagonal elements and multiply non-diagonal elements by dt
    for i in 1:size(k_matrix, 1)
        for j in 1:size(k_matrix, 2)
            if i == j
                k_matrix[i, i] = 1 - dt * sum(k[i, [1:i-1; i+1:end]])
            else
                k_matrix[i, j] *= dt
            end
        end
    end
    return k_matrix
end

function arrhenius_transform_with_deposition(existing_layers, new_layer, k_m) #in the time step of dt
    # Transform existing layers
    transformed_layers = existing_layers * k_m
    transformed_layers = vcat(transformed_layers, reshape(new_layer, 1, :))

    return transformed_layers
end