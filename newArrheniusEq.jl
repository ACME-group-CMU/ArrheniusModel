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

function flow_coefficient(type, effecting_nums, decay_coefficient)
"""
Inputs:
    type: type of decay (exponential or linear)
    effecting_nums: number of effecting layers on top of the current layer
Output:
    flow_coefficients: array of flow coefficients for each layers
"""
    effecting_nums = round(Int, effecting_nums)
    fcoeff = zeros(effecting_nums)
    for j in 1:effecting_nums
        if type == "exponential"
            fcoeff[j] = exp(-decay_coefficient * j)
        elseif type == "linear"
            fcoeff[j] = (1 - decay_coefficient * j)
        end
    end
    return fcoeff
end

function arrhenius_rate(pe::PhaseEnergies, T=300)
    kb = 8.617e-5 #eV/K
    A = 1.0 # Arrhenius prefactor
    pe.K = A * exp.(-pe.Ea_plus_ΔG ./ (kb * T))
    # Adjust the diagonal elements
    for i in 1:size(pe.K, 1)
        pe.K[i, i] =  -1 * sum(pe.K[i, [1:i-1; i+1:end]])
    end
end