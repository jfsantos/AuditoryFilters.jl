using Auditory, DSP, Base.Test

s = read(open(joinpath(dirname(@__FILE__), "data", "sa2.raw")), Float64, 42701)
fs = 16000
(W, _) = Auditory.fft2gammatonemx(512, fs, 23, 1, 150, 4000)
G, F = gammatonegram(s, fs, 0.025, 0.010, 23, 150, 4000, 1)

G_matlab = readcsv(open(joinpath(dirname(@__FILE__), "data", "G_matlab.csv")))
F_matlab = readcsv(open(joinpath(dirname(@__FILE__), "data", "F_matlab.csv")))
W_matlab = readcsv(open(joinpath(dirname(@__FILE__), "data", "W_matlab.csv")))

@test_approx_eq W W_matlab
@test_approx_eq F F_matlab'
@test_approx_eq G G_matlab
