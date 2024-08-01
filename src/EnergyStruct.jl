using StaticArrays
struct PhaseEnergies
    G::AbstractVector
    Ea::AbstractMatrix
    barriers::AbstractMatrix
    function PhaseEnergies(G::AbstractVector, forward_Ea::AbstractMatrix)
        n = length(G)
        @assert size(forward_Ea) == (n, n)
        deltaG = ΔG(G)
        barriers = [i == j ? 0 : (deltaG[i, j] > 0 ? deltaG[i, j] + forward_Ea[i, j] : forward_Ea[i, j])
                for j in 1:n, i in 1:n]
        new(G, forward_Ea, Matrix{eltype(G)}(barriers))
    end
end

n_phases(pe::PhaseEnergies) = length(pe.G)
ΔG(gmat::AbstractVector) = [gmat[j] - gmat[i] for j in eachindex(gmat), i in eachindex(gmat)]
ΔG(pe::PhaseEnergies) = ΔG(pe.G)