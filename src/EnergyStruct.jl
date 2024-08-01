using StaticArrays

mutable struct PhaseEnergies{N,F}
    G::SVector{N,F}
    Ea::SMatrix{N,N,F}
    ΔG::SMatrix{N,N,F}
    Ea_plus_ΔG::SMatrix{N,N,F}
    function PhaseEnergies(G::AbstractVector, Ea::AbstractMatrix)
        n = length(G)
        F = eltype(G)
        @assert size(Ea) == (n, n)
        ΔG = [G[j] - G[i] for j in 1:n, i in 1:n]
        Ea_plus_ΔG = [i == j ? 0 : (ΔG[i, j] > 0 ? ΔG[i, j] + Ea[i, j] : Ea[i, j])
                for j in 1:n, i in 1:n]
        new{n,F}(SVector{n,F}(G), SMatrix{n,n,F}(Ea), SMatrix{n,n,F}(ΔG), SMatrix{n,n,F}(Ea_plus_ΔG))
    end
end

n_phases(pe::PhaseEnergies) = length(pe.G)
