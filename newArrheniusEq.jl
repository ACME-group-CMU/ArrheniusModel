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
"""
function flow_coefficient(type, effecting_layers, decay_constant)
    A = 1.0
    if type == "exponential"
        return A * exp(-decay_constant * coveringlayers)
    elseif type == "linear"
        return A * (1 - decay_constant * coveringlayers)
    end
end
"""

function flow_coefficient(type, effecting_nums, decay_coefficient)
    # effecting_nums: it approaches zero and stop decreasing as the number of effecting layers increases to a certain point, 
    #                 so I manually set this to limit the number of effecting layers
    flow_coefficients = []
    for i in 1:effecting_nums
        if type == "exponential"
            push!(flow_coefficients, exp(-decay_coefficient * i))
        elseif type == "linear"
            push!(flow_coefficients, (1 - decay_coefficient * i))
        end
    end
    return flow_coefficients
end

function arrhenius_rate(Ea, T)
    k = exp(-Ea / (kb* T))
    return k
end

function arrhenius_rate_matrix(Ea_m, T)
    k_matrix = arrhenius_rate.(Ea_m, T)
    k = copy(k_matrix)
    # Adjust the diagonal elements
    for i in 1:size(k_matrix, 1)
        k_matrix[i, i] =  -1 * sum(k[i, [1:i-1; i+1:end]])
    end
    return k_matrix
end
