
using Test, SHA, SparseArrays
using SBML
import Pkg

@testset "SBML test suite" begin
    include("version.jl")
    include("ecoli_flux.jl")
    include("ecoli_modelonly.jl")
end
