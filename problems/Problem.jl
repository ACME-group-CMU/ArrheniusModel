using DifferentialEquations, Enzyme, ArrheniusModel

G = [-5.92, -5.942, -5.97]
Ea = [0.00 1.00 0.01; 1.00 0.00 1.00; 0.01 1.00 0.00]
pe = PhaseEnergies(G, Ea)
T = 300.0
t= 10
dt = 0.05
datasize = Int(t/0.5+1)
num_steps = floor(Int, t/dt)
num_layers = floor(Int, t/0.5)+1
flow_rate = 0.5
decay_coefficient = 0.00001 * flow_rate
fcoeff = flow_coefficient("exponential", num_layers, decay_coefficient)
n = n_phases(pe)
c0 = zeros(num_layers, n)
c0[1, 1] = 1.0
j = 0
j0 = 0
p = (fcoeff, j0, j, dt, num_steps, num_layers)
tspan = (0.0, (num_steps-1) * dt)
tsteps = range(tspan[1], tspan[2]; length = datasize)

function deposition_rates!(dc, c, p, t)
    # Unpack parameters
    fcoeff, j0, j, dt, num_steps, num_layers = p
    # Calculate deposition rates
    j = floor(Int, t / 0.5) + 1
    f = reverse(fcoeff[j: num_layers+j-1])
    dc .= c .* f * K
    if j != j0
        c[j+1, 1] = 1.0
        j = j0
    end
end

function ODE_calculation(barriers=pe.barriers, T=T) #defining the prob within the function
    K = arrhenius_rate(barriers, T)
    prob = ODEProblem(deposition_rates!, c0, tspan, p)
    ode_data = Array(solve(prob, Euler(), saveat = 0.5, dt = 0.05))
    return ode_data[:,:,end]
end

db = Array(zero(pe.barriers))
dT = 1.0

godeT = Enzyme.autodiff(Forward,ODE_calculation,Duplicated(pe.barriers, db), Duplicated(T,dT))[1]
#The godeT shows the ODE_data instead of its derivative wrt T
display(ODE_calculation())
display(godeT)


prob = ODEProblem(deposition_rates!, c0, tspan, p) #defining the prob out of the function
function ODE_calculation(barriers=pe.barriers, T=T, prob=prob)
    K = arrhenius_rate(barriers, T)
    ode_data = Array(solve(prob, Euler(), saveat = 0.5, dt = 0.05))
    return ode_data[:,:,end]
end

db = Array(zero(pe.barriers))
dT = 1.0

godeT = Enzyme.autodiff(Forward,ODE_calculation,Duplicated(pe.barriers, db), Duplicated(T,dT), Const(prob))[1]
#ERROR: LoadError: Constant memory is stored (or returned) to a differentiable variable.
display(ODE_calculation())
display(godeT)