module ArrheniusModel

using Statistics
using LinearAlgebra
using DifferentialEquations

include("EnergyStruct.jl")
include("newArrheniusEq.jl")
include("Energies.jl")
include("Simulation.jl")

export PhaseEnergies, n_phases
export flow_coefficient, meshgrid
export deposition_rates!, simulate_deposition, most_preferable_state, arrhenius_rate
# TODO: add exports for anything else you would want to use extnerally, i.e. not functions that are just internal "helper" things

end
