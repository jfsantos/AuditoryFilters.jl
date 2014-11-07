using Auditory, DSP, Base.Test

s = read(open(joinpath(dirname(@__FILE__), "data", "sa2.raw")), Float64, 42701)
fs = 16000
(W, _) = Auditory.fft2gammatonemx(512, fs, 23, 1, 150, 4000)
nfft = 512
nwin = int(0.025*fs)
nhop = int(0.010*fs)
S = stft(s, nwin, nwin-nhop; nfft=nfft, fs=fs, window=hanning)
G, F = gammatonegram(s, fs, 0.025, 0.010, 23, 150, 4000, 1)

S_matlab_re = readcsv(open(joinpath(dirname(@__FILE__), "data", "S_matlab_re.csv")))
S_matlab_im = readcsv(open(joinpath(dirname(@__FILE__), "data", "S_matlab_im.csv")))
S_matlab = complex(S_matlab_re, S_matlab_im)
G_matlab = readcsv(open(joinpath(dirname(@__FILE__), "data", "G_matlab.csv")))
F_matlab = readcsv(open(joinpath(dirname(@__FILE__), "data", "F_matlab.csv")))
W_matlab = readcsv(open(joinpath(dirname(@__FILE__), "data", "W_matlab.csv")))

@test_approx_eq W W_matlab
@test_approx_eq F F_matlab'
@test_approx_eq G G_matlab
