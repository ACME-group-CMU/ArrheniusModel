using SciMLSensitivity, OrdinaryDiffEq, Enzyme

function fiip(du, u, p, t)
    du[1] = dx = p[1] * u[1] - p[2] * u[1] * u[2]
    du[2] = dy = -p[3] * u[2] + p[4] * u[1] * u[2]
end
p = [1.5, 1.0, 3.0, 1.0];
u0 = [1.0; 1.0];
prob = ODEProblem(fiip, u0, (0.0, 10.0), p)
sol = solve(prob, Tsit5())[end]
print(sol)
#loss(u0, T, p, prob) = sum(solve(prob, Tsit5(), u0 = u0, p = p, saveat = 0.1))
last(u0, p, prob) = solve(prob, Tsit5(), u0 = u0, p = p, saveat = 0.1)[end]
function last!(u0, T, p, prob)
    u0 .= last(u0, p, prob)
    return nothing
end
a = last!(u0, T, p, prob)
println(a, u0)
du = zeros(size(u0))
dp = zeros(size(p))
du[1] = 1.0
T = [1.0]
dT = zero(T)
Enzyme.autodiff(ReverseWithPrimal, last!, Duplicated(u0, du), Duplicated(T, dT), Duplicated(p, dp), Const(prob))