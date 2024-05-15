using ArrheniusModel
using Test

@testset "ArrheniusModel.jl" begin
    # Write your tests here.
    # divide tests up into multiple subsets, @testset's can be nested
    # here's an example, feel free to add to/change it
    @testset "PhaseEnergies struct" begin
        G = [0.0,0.0]
        Ea = [0. 0.2; 0.2 0.]
        pe = PhaseEnergies(G, Ea)
        @test n_phases(pe) == 2
        @test pe.Ea_plus_ΔG == [0.0 0.2; 0.2 0.0]
        @test_throws AssertionError PhaseEnergies([0,0,0],Ea)
    end
    @testset "Homogeneous 3" begin
        G = [0.0,0.0,0.0]
        Ea = [0. 0. 0.; 0. 0. 0.; 0. 0. 0.]
        pe = PhaseEnergies(G, Ea)
        @test n_phases(pe) == 3
        @test all(pe.Ea_plus_ΔG .== 0)
        @testset "300K simulation" begin
            T = 300.0
            num_steps = 10
            dt = 0.1
            flow_rate = 0.5
            decay_coefficient = 0.00001 * flow_rate
            fcoeff = flow_coefficient("exponential", num_steps, decay_coefficient)
            layers = simulate_deposition(fcoeff, pe, T, num_steps, dt)
            @test layers[:, 2] == layers[:, 3]
            @test size(layers) == (num_steps+1, 3)
            @test all(layers .>= 0)
            @test all(layers[1,:] .== [1.0, 0.0, 0.0])
            phase = most_preferable_state(layers, 0.01, ["A", "B", "C"])
            @test phase == "A+B+C"
            phase = most_preferable_state(layers, 0.3, ["A", "B", "C"])
            @test phase == "A"
        end
    end
    @testset "Low K Test" begin
        G = [-5.92, -5.942, -5.97]
        Ea = [0.00 1.00 0.01; 1.00 0.00 1.00; 0.01 1.00 0.00]
        pe = PhaseEnergies(G, Ea)
        T = 1.0
        num_steps = 10
        dt = 0.1
        flow_rate = 0.5
        decay_coefficient = 0.00001 * flow_rate
        fcoeff = flow_coefficient("exponential", num_steps, decay_coefficient)
        layers = simulate_deposition(fcoeff, pe, T, num_steps, dt)
        @test size(layers) == (num_steps+1, 3)
        @test all(layers .>= 0)
        @test all(layers[1,:] .== [1.0, 0.0, 0.0])
        @test all(layers[:, 1] .== 1.0)
        for i in 1:size(pe.K,1)
            @test pe.K[i, i] ==  -1 * sum(pe.K[i, [1:i-1; i+1:end]])
        end
    end
    @testset "this will not fail" begin
        @test 1==1
    end
end
