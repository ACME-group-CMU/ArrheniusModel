using StaticArrays

mutable struct PhaseEnergies{N,T}
    G::SVector{N,T}
    Ea::SMatrix{N,N,T}
    ΔG::SMatrix{N,N,T}
    Ea_plus_ΔG::SMatrix{N,N,T}
    K::MMatrix{N,N,T}
    function PhaseEnergies(G::AbstractVector, Ea::AbstractMatrix)
        n = length(G)
        T = eltype(G)
        @assert size(Ea) == (n, n)
        ΔG = [G[j] - G[i] for j in 1:n, i in 1:n]
        Ea_plus_ΔG = [i == j ? 0 : (ΔG[i, j] > 0 ? ΔG[i, j] + Ea[i, j] : Ea[i, j])
                for j in 1:n, i in 1:n]
        new{n,T}(SVector{n,T}(G), SMatrix{n,n,T}(Ea), SMatrix{n,n,T}(ΔG), SMatrix{n,n,T}(Ea_plus_ΔG), MMatrix{n,n,T}(zeros(n,n)))
    end
end

n_phases(pe::PhaseEnergies) = length(pe.G)
