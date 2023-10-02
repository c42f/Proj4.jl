using Proj
using Test
using StaticArrays: SA
import Aqua

# Base.isapprox doesn't work on tuples
is_approx(a, b) = all(isapprox(c[1], c[2]) for c in zip(a, b))

@testset "Proj" begin
    include("libproj.jl")
    include("applications.jl")
    Aqua.test_all(Proj; ambiguities = false)
end # testset "Proj"
