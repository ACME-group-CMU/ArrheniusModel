using SciMLSensitivity, OrdinaryDiffEq, Enzyme
using ArrheniusModel

T = 573.0
flow_rate = 0.5

#change the barriers to 3*3
barriers = [0.0 0.2 0.4; 0.2 0.0 0.3; 0.4 0.3 0.0]
#barriers = [0.0 0.2; 0.2 0.0]
K = arrhenius_rate(barriers, T)
#fcoeff = flow_coefficient.("exponential", num_layers, decay_constant*flow_rate)
n = size(barriers, 1)
t = 60 # seconds
dt = 0.5 # seconds
num_steps = round(Int, t/dt)
num_layers = floor(Int, t/0.5)+1
j = 0
j0 = 0
tspan = (0.0, (num_steps-1) * dt)
decay_constant = 0.01

function fiip(dc, c, p, t)
    #T, flow_rate = p
    K = arrhenius_rate(barriers, p[1])
    fcoeff = flow_coefficient("exponential", num_layers, decay_constant*p[2])
    j = floor(Int, t / 0.5) + 1
    f = reverse(fcoeff[j: num_layers+j-1])
    #f = ones(num_layers,1) .* exp.(-(1:num_layers) * p[2])
    """
    When I use dc .= c .*f * K, it gives me activity error
    Constant memory is stored (or returned) to a differentiable variable.
    """
<<<<<<< HEAD
    dc .= c * K .* f
    #dc .= c .*f * K
=======
    # dc .= c * K .*f
    dc .= c .*f * K
>>>>>>> 27c41f595e243d6b7154b2a438519410f8d34829
    if j != j0
        c[j+1, 1] = 1.0
        j = j0
    end
end
p = [T, flow_rate]
c0 = zeros(num_layers, n)
c0[1,1] = 1.0
prob = ODEProblem(fiip, c0, tspan, p)
#sol = solve(prob, Euler(), dt=0.1)[end]
#display(sol)
#loss(u0, T, p, prob) = sum(solve(prob, Tsit5(), u0 = u0, p = p, saveat = 0.1))
last(u0, p, prob) = solve(prob, Euler(), u0 = u0, p = p, saveat = 0.1, dt=0.1)[end]
display(last(c0, p, prob))

function last!(u0, p, prob)
    u0 .= last(u0, p, prob)
    return nothing
end
"""
function last!(u0, T, f, prob)
    u0 .= last(u0, [T[1], f[1]], prob)
    return nothing
end
"""
#a = last!(c0, T, flow_rate, prob)
#display(a)
#display(c0)
dc0 = zeros(size(c0))
dc0[2,2] = 1.0
dp = zeros(size(p))
db = zeros(size(barriers))
f_ = [flow_rate]
df = zero(f_)
T_ = [T]
dT = zero(T_)
#Enzyme.autodiff(Reverse, last!, Duplicated(c0, dc0), Duplicated(T_, dT), Duplicated(f_, df), Const(prob))
Enzyme.autodiff(Reverse, last!, Duplicated(c0, dc0), Duplicated(p,dp), Const(prob))



