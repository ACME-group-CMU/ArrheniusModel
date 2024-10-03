using SciMLSensitivity, OrdinaryDiffEq, Enzyme
using ArrheniusModel
"""
1. I don't know why i have to assign a value to du to obtain the gradient, what does the du mean?
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
#display(last(u0, p, prob))
du = zeros(size(u0))
dp = zeros(size(p))
du[1]=1.0
#println(sol)
gloss = Enzyme.autodiff(Reverse, last, Duplicated(u0, du), Duplicated(p, dp), Const(prob))
println("Automatic gradient", dp)

# another attempt, let's make it single vector input/scalar output for simplicity...
function last3(u0, p, prob)
    solve(prob, Tsit5(), u0 = [u0], p = [p], saveat = 0.1)[end][1]
end

# first, what about the convenience function?
gradient(Reverse, x -> last3(x[1],x[2],prob), [1.0, 1.0])
# okay, so this seems to work? How about the other way...

u0 = 1.0
p = 1.0
Enzyme.autodiff(Reverse, last3, Active(u0), Active(p), Const(prob))
# also good...okay, interesting...

# here's a thought...
function last4(u0, p, prob, val)
    val .= solve(prob, Tsit5(), u0 = u0, p = p, saveat = 0.1)[end]
    return nothing
end
p = [1.0]
u0 = [1.0]
du = zeros(size(u0))
dp = zeros(size(p))
val = [0.0]
dval = zeros(size(val))
dval[1] = 1.0
# okay, I think this is now correct? But the question is still, what exactly was wrong before?
# that is, why does overwriting u0 in that way cause the gradient to square? Mysterious scope things?

#Finite_difference to verify
function last2(u0, p, prob) #returning the last value directly
    return solve(prob, Tsit5(), u0 = u0, p = p, saveat = 0.1)[end]
end
Enzyme.autodiff(Reverse, last4, Duplicated(u0, du), Duplicated(p, dp), Const(prob), Duplicated(val, dval))

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


"""
2. The diagonal value is missing in Reverse mode, why? Also, dT doesn't have any effect on the result in this form
"""
# example from tests
G = [1.0,0.0]
Ea = [0. 0.2; 0.2 0.]
pe = PhaseEnergies(G, Ea)
#Forward method
db = Array(zero(pe.barriers))
db[1,2] = 1.0
g12 = Enzyme.autodiff(Forward, arrhenius_rate, Duplicated(pe.barriers, db))[1]
println(g12)
#Reverse method
db = zero(pe.barriers)
dT = 0.0
T = 300.0
K = zero(pe.barriers)
dK = one(pe.barriers)
# dK[1,2] = 1.0
dK[1,1] = 1.0
Enzyme.autodiff(ReverseWithPrimal, arrhenius_rate!, Duplicated(pe.barriers, db), Duplicated(K, dK), Duplicated(T, dT))
println(db) #The diagonal value is missing, why? Also, dT doesn't have any effect on the result

# first, a version that takes in a matrix and returns one of the same size
function dummy_ar2!(b,T,val)
    val .= [arrhenius_rate(b, T[1])[1,1]]
    return nothing
end
T = [300.0]
dT = [0.0]
db = zero(pe.barriers)
val = [0.0]
dval = [1.0]
Enzyme.autodiff(ReverseWithPrimal, dummy_ar2!, Duplicated(pe.barriers, db), Duplicated(T, dT), Duplicated(val,dval))

# now let's have it return the full matrix...
ar_matrixT!(b,K,T) = arrhenius_rate!(b, K, T[1])
T = [300.0]
dT = [0.0]
# b = [0.0 1.0; 1.2 0.0]
b = [0 0.5; 0.5 0]
db = zero(b)
K = zero(b)
dK = zero(b)
dK[1,2] = 1.0
# dK[2,1] = 1.0
Enzyme.autodiff(ReverseWithPrimal, ar_matrixT!, Duplicated(b, db), Duplicated(K,dK), Duplicated(T,dT))
# so we don't have to change the actual package functions, just whenever we want T-sensitivity, we have to make sure to pass it into Enzyme as a 1 x 1 matrix