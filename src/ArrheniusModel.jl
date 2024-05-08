module ArrheniusModel

include("EnergyStruct.jl")
include("simulation.jl")
# TODO: add the rest of them

export PhaseEnergies, n_phases
export simulate_deposition, most_preferable_state
# TODO: add exports for anything else you would want to use extnerally, i.e. not functions that are just internal "helper" things

end
