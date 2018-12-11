# Testing modulation filterbank response
using Compat.DelimitedFiles
x = [1, zeros(9999)...]
mf = compute_modulation_cfs(4,128,8)
mfb = modulation_filterbank(mf, 16000, 2)
y = filt(mfb, x)

mf_matlab = readdlm(open(joinpath(dirname(@__FILE__), "data", "modfilters_freqs.csv")), ',')
y_matlab = readdlm(open(joinpath(dirname(@__FILE__), "data", "modfilters_response.csv")), ',')

@test y ≈ y_matlab'
@test mf ≈ mf_matlab
