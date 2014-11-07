using Auditory, Base.Test

x = [1, zeros(9999)]
mf = compute_modulation_cfs(4,128,8)
mfb = modulation_filterbank(mf, 16000, 2)
y = filt(mfb, x)

mf_matlab = readcsv(open(joinpath(dirname(@__FILE__), "data", "modfilters_freqs.csv")))
y_matlab = readcsv(open(joinpath(dirname(@__FILE__), "data", "modfilters_response.csv")))

@test_approx_eq y y_matlab'
@test_approx_eq mf mf_matlab
