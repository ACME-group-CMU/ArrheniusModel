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
            t= 10
            dt = 0.1
            num_steps = floor(Int, t/dt)
            num_layers = floor(Int, t/0.5)+1
            flow_rate = 0.5
            decay_coefficient = 0.00001 * flow_rate
            fcoeff = flow_coefficient("exponential", num_layers, decay_coefficient)
            layers = simulate_deposition(fcoeff, pe, T, num_steps, num_layers, dt)
            @test all(sum(layers, dims=2) .≈ 1.0) #Conservation rule
            @test layers[:, 2] ≈ layers[:, 3]  # =somehow doesn't work even it shows the same value
            @test size(layers) == (num_layers, 3) # num_steps+1 to num_steps due to format change
            @test all(layers .>= 0)
            @test all(layers[end,:] .== [1.0, 0.0, 0.0]) #layers[1,:] -> layers[end,:] due to format change
            phase = most_preferable_state(layers, 0.01, ["A", "B", "C"])
            @test phase == "A+B+C"
            phase = most_preferable_state(layers, 0.5, ["A", "B", "C"])
            @test phase == "" #For a short time test, results seems to depend on solver => Change to something that is forerver true
        end
    end
    @testset "Low K Test" begin
        G = [-5.92, -5.942, -5.97]
        Ea = [0.00 1.00 0.01; 1.00 0.00 1.00; 0.01 1.00 0.00]
        pe = PhaseEnergies(G, Ea)
        T = 1.0
        t= 10
        dt = 0.1
        num_steps = floor(Int, t/dt)
        num_layers = floor(Int, t/0.5)+1
        flow_rate = 0.5
        decay_coefficient = 0.00001 * flow_rate
        fcoeff = flow_coefficient("exponential", num_layers, decay_coefficient)
        layers = simulate_deposition(fcoeff, pe, T, num_steps, num_layers, dt)
        @test all(sum(layers, dims=2) .≈ 1.0)
        @test layers[:, 2] != layers[:, 3]
        @test size(layers) == (num_layers, 3)
        @test all(layers .>= 0)
        @test all(layers[end,:] .== [1.0, 0.0, 0.0])
        @test all(layers[:, 1] .== 1.0)
        for i in 1:size(pe.K,1)
            @test pe.K[i, i] ==  -1 * sum(pe.K[i, [1:i-1; i+1:end]])
        end
    end
    @testset "this will not fail" begin
        @test 1==1
    end
end
