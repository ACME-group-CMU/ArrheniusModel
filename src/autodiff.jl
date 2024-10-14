using ChainRules
using Enzyme
using ArrheniusModel

# example from tests
G = [1.0,0.0]
Ea = [0. 0.2; 0.2 0.]
pe = PhaseEnergies(G, Ea)

# try forward mode without the struct first
# this will give the derivative of each element of the output matrix wrt the 1,2 element of the input
db = Array(zero(pe.barriers))
db[1,2] = 1.0
g12 = Enzyme.autodiff(Forward, arrhenius_rate, Duplicated(pe.barriers, db))[1]
# (convince yourself that this gives the right result)

# similarly, to get the derivative wrt the 2,2 element we can do...
db = Array(zero(pe.barriers))
db[2,2] = 1.0
Enzyme.autodiff(Forward, arrhenius_rate, Duplicated(pe.barriers, db))
# (again, convince yourself)

# we can do it wrt the temperature also...
db = Array(zero(pe.barriers))
dT = 1.0
Enzyme.autodiff(Forward, arrhenius_rate, Duplicated(pe.barriers, db), Duplicated(300.0, dT))
# (check this on paper too)

# let's try reverse mode...
db = Array(zero(pe.barriers))
db[1,2] = 1.0
Enzyme.autodiff(Reverse, arrhenius_rate, Duplicated(pe.barriers, db))
# this doesn't work "as-is" because for reasons I don't entirely understand (but are documented here: https://enzyme.mit.edu/julia/stable/api/#EnzymeCore.autodiff-Union{Tuple{Nargs},%20Tuple{Holomorphic},%20Tuple{RABI},%20Tuple{ReturnPrimal},%20Tuple{A},%20Tuple{FA},%20Tuple{ReverseMode{ReturnPrimal,%20RABI,%20Holomorphic},%20FA,%20Type{A},%20Vararg{Annotation,%20Nargs}}}%20where%20{FA%3C:Annotation,%20A%3C:Annotation,%20ReturnPrimal,%20RABI%3C:EnzymeCore.ABI,%20Holomorphic,%20Nargs}), you can't "directly" use this with functions that don't return scalar values, but instead have to turn them into mutating functions...so let's try that
db = zero(pe.barriers)
K = zero(pe.barriers)
dK = copy(K)
dK[1,2] = 1.0
Enzyme.autodiff(Reverse, arrhenius_rate!, Duplicated(pe.barriers, db), Duplicated(K, dK))
# the difference now is that our "result" shows up in `db` and it's the `db` value that would have led to the `dK` value we set (or something proportional to that anyway...I'm still a little fuzzy)...we can check that this is all consistent by feeding in the output gradient from the forward-mode version and making sure we just get a 1 in `db`...
db = zero(pe.barriers)
K = zero(pe.barriers)
Enzyme.autodiff(Reverse, arrhenius_rate!, Duplicated(pe.barriers, db), Duplicated(K, g12))
# okay so we didn't get a 1 there, but we did get something proportional to 1 so maybe there's an additional chain rule to account for that I'm not thinking through clearly just yet, so some more pen-and-paper checks are definitely needed

# to differentiate "wrt" our struct, we need to get a bit fancier and define a custom rule...
# I HAVEN'T FINISHED WITH THIS YET AND IT DEFINITELY WON'T WORK RIGHT NOW
function ChainRulesCore.rrule(::typeof(arrhenius_rate), pe::PhaseEnergies, T::Real) 
    K = arrhenius_rate(pe, T)
    function arrhenius_rate_pullback(K̄)
        f̄ = NoTangent()
        kb = 8.617e-5
        p̄e = Tangent{PhaseEnergies}(; G=NoTangent(), Ea=NoTangent(), barriers=-(K̄./K̄)/(kb*T))
        T̄ = @thunk(foo.A' * ȳ)
        return f̄, p̄e, T̄
    end
    return K, foo_mul_pullback
end