using ArrheniusModel
using Test

@testset "ArrheniusModel.jl" begin
    # Write your tests here.
    # divide tests up into multiple subsets, @testset's can be nested
    # here's an example, feel free to add to/change it
    @testset "PhaseEnergies struct" begin
        G = [0.0,0.0]
        Ea = [0.1 0.2; 0.3 0.4]
        pe = PhaseEnergies(G, Ea)
        @test n_phases(pe)==2
        @test_throws AssertionError PhaseEnergies([0,0,0],Ea)
    end
end
