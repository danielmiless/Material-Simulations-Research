using Pkg
Pkg.activate(@__DIR__ * "/..")

using Test

# Test suite for Material Simulations
# Add tests here as the project develops

@testset "Material Simulations Tests" begin
    @test 1 + 1 == 2
    
    # Add more tests here
    # Example:
    # @testset "Exponential Spring Model" begin
    #     # Test exponential spring calculations
    # end
end

