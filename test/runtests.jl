using AuditoryFilters, DSP, Compat, Compat.Test

@testset "ERB Filterbank" begin
    include("test_erb_fb.jl")
end

@testset "Modulation Filterbank" begin
    include("test_modulation_fb.jl")
end

@testset "Gammatonegram" begin
    include("test_gammatonegram.jl")
end
