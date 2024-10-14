using ArrheniusModel
using Statistics
using DifferentialEquations
using Plots
using Printf
using Enzyme
gr()
pe = PhaseEnergies(ArrheniusModel.G_values, ArrheniusModel.Ea_constants)
display(pe.barriers)

#Run the single point Simulation
T = 573.0
flow_rate = 0.5
threshold = 0.3 # Threshold for most preferable state
t = 60 # seconds
dt = 0.5 # seconds
num_steps = round(Int, t/dt)
num_layers = floor(Int, t/0.5)+1
para_sim = num_steps, num_layers, dt
phase_names = ["x", "α", "β", "κ", "γ", "δ"]
compositions_all = simulate_deposition(flow_rate, T, pe.barriers, para_sim)
display(compositions_all)

#SensitivityAnalysis FAILED!!
sens_f, sens_T, sens_b = sens_simulation(flow_rate, T, pe, para_sim)
display(sens_f)
display(sens_T)
display(sens_b)