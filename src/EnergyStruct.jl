using StaticArrays
struct PhaseEnergies{N,F}
    G::SVector{N,F}
    Ea::SMatrix{N,N,F}
    barriers::SMatrix{N,N,F}
    function PhaseEnergies(G::AbstractVector, forward_Ea::AbstractMatrix)
        n = length(G)
        F = eltype(G)
        @assert size(forward_Ea) == (n, n)
        deltaG = ΔG(G)
        barriers = [i == j ? 0 : (deltaG[i, j] > 0 ? deltaG[i, j] + forward_Ea[i, j] : forward_Ea[i, j])
                for j in 1:n, i in 1:n]
        new{n,F}(SVector{n,F}(G), SMatrix{n,n,F}(forward_Ea), SMatrix{n,n,F}(barriers))
    end
end

n_phases(pe::PhaseEnergies{N,F}) where {N,F} = N
ΔG(gmat::AbstractVector) = [gmat[j] - gmat[i] for j in eachindex(gmat), i in eachindex(gmat)]
ΔG(pe::PhaseEnergies) = ΔG(pe.G)