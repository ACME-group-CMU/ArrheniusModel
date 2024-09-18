using SciMLSensitivity, OrdinaryDiffEq, Enzyme
"""
1. I don't know why i have to assign a value to du to obtain the gradient, what does the du mean?
2. The gradient from AD is the square value according to the FD method, why?
Symbolic method for this case:
du/dt = u*p
u^-1 du = p dt
ln(u) = p*t
u = exp(p*t)
du/dp = t*exp(p*t) 
In this case, t = 1.0, p = 1.0, du/dp = e, consistent to FD but not AD
"""


#AD test
function simpletest(du, u, p, t)
    du[1] = u[1] * p[1]
end
p = [1.0]
u0 = [1.0]
prob = ODEProblem(simpletest, u0, (0.0, 1.0), p)
sol = solve(prob, Tsit5(), saveat =0.1)
#display(sol)
function last(u0, p, prob)
    u0 .= solve(prob, Tsit5(), u0 = u0, p = p, saveat = 0.1)[end]
    return nothing
end
display(last(u0, p, prob))
du = zeros(size(u0))
dp = zeros(size(p))
du[1]=1.0
#println(sol)
gloss = Enzyme.autodiff(Reverse, last, Duplicated(u0, du), Duplicated(p, dp), Const(prob))
println(dp)


#Finite_difference to verify
function last2(u0, p, prob) #returning the last value directly
    return solve(prob, Tsit5(), u0 = u0, p = p, saveat = 0.1)[end]
end

function finite_difference_gradient(f, u0, p, prob, epsilon=1e-6)
    grad = zeros(length(p))
    for i in 1:length(p)
        p_forward = copy(p)
        p_backward = copy(p)
        p_forward[i] += epsilon
        p_backward[i] -= epsilon
        a = f(u0, p_forward, prob)
        b = f(u0, p_backward, prob)
        grad[i] = (a[1]-b[1]) / (2 * epsilon)
    end
    return grad
end

# Reset u0 to its initial value
u0 = [1.0]

# Compute the numerical gradient
numerical_grad = finite_difference_gradient(last2, u0, p, prob)
println("Numerical Gradient: ", numerical_grad)