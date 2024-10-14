module ArrheniusModel

using Statistics
using LinearAlgebra
using DifferentialEquations

include("EnergyStruct.jl")
include("newArrheniusEq.jl")
include("Energies.jl")
include("Simulation.jl")
include("SensitivityAnalysis.jl")

export PhaseEnergies, n_phases
export flow_coefficient, meshgrid, arrhenius_rate, arrhenius_rate!, ar_matrixT!
export deposition_rates!, simulate_deposition, most_preferable_state, simulate_deposition!
export sens_arrhenius_rate, sens_simulation
# TODO: add exports for anything else you would want to use extnerally, i.e. not functions that are just internal "helper" things

end
