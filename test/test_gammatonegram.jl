# Testing gammatonegram response
using Compat.DelimitedFiles
s = read!(open(joinpath(dirname(@__FILE__), "data", "sa2.raw")), Array{Float64}(undef,42701))
fs = 16000
(W, _) = AuditoryFilters.fft2gammatonemx(512, fs, 23, 1, 150, 4000)
G = gammatonegram(s, fs, 0.025, 0.010, 23, 150, 4000, 1)

G_matlab = readdlm(open(joinpath(dirname(@__FILE__), "data", "G_matlab.csv")), ',')
F_matlab = readdlm(open(joinpath(dirname(@__FILE__), "data", "F_matlab.csv")), ',')
W_matlab = readdlm(open(joinpath(dirname(@__FILE__), "data", "W_matlab.csv")), ',')

@test W ≈ W_matlab
@test G.frequencies ≈ F_matlab
@test G.amplitude ≈ G_matlab
