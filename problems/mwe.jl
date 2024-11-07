using SciMLSensitivity, OrdinaryDiffEq, Enzyme

# this one works, AD runs, but only gives sensitivity for p[1]
function deriv_v1(dy, y, p, t)
    M = [3 1 4; 1 5 9; 6 5 4] .^ p[1]
    f = ones(4,1) .* exp.(-[1, 2, 3, 4] * p[2])
    dy .= y * M .* f
end

# this one works but AD gives activity error
function deriv_v2(dy, y, p, t)
    M = [3 1 4; 1 5 9; 6 5 4] .^ p[1]
    f = ones(4,1) .* exp.(-[1, 2, 3, 4] * p[2])
    dy .= y .* f * M # this line is different
end

p = [0.42, 7]
y0 = ones(4,3)
tspan = (0.0, 1.0)
prob1 = ODEProblem(deriv_v1, y0, tspan, p)
prob2 = ODEProblem(deriv_v2, y0, tspan, p)

last(u0, p, prob) = solve(prob, Euler(), u0 = u0, p = p, saveat=0.1, dt=0.1)[end]
display(last(y0, p, prob1))
display(last(y0, p, prob2))

function last!(u0, p, prob)
    u0 .= last(u0, p, prob)
    return nothing
end

dy01 = zeros(size(y0))
dy01[2,2] = 1.0
dp1 = zeros(size(p))
Enzyme.autodiff(Reverse, last!, Duplicated(y0, dy01), Duplicated(p,dp1), Const(prob1))

dy02 = zeros(size(y0))
dy02[2,2] = 1.0
dp2 = zeros(size(p))
Enzyme.autodiff(Reverse, last!, Duplicated(y0, dy02), Duplicated(p,dp2), Const(prob2))