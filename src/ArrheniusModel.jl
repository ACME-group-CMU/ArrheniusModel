module ArrheniusModel

using Statistics
using LinearAlgebra

include("EnergyStruct.jl")
include("newArrheniusEq.jl")
include("Energies.jl")
include("Simulation.jl")

export PhaseEnergies, n_phases
export simulate_deposition, most_preferable_state, flow_coefficient, meshgrid
# TODO: add exports for anything else you would want to use extnerally, i.e. not functions that are just internal "helper" things

end
