using Auditory, DSP, MATLAB, Base.Test, WAV

(s, fs) = wavread("/Users/jfsantos/Projects/SRMR_toolbox/gammatonegram/sa2.wav")
s = s[:,1]
(W, _) = Auditory.fft2gammatonemx(512, fs, 23, 1, 150, 4000)
nfft = 512
nwin = int(0.025*fs)
nhop = int(0.010*fs)
S = stft(s, nwin, nwin-nhop; nfft=nfft, fs=fs, window=hanning)
G, F = gammatonegram(s, fs, 0.025, 0.010, 23, 150, 4000, 1)

fs = float(fs)
nfft = float(nfft)
nwin = float(nwin)
nhop = float(nhop)
restart_default_msession()
@mput s
@mput fs
@mput nfft
@mput nwin
@mput nhop
@matlab begin
	addpath("/Users/jfsantos/Projects/SRMR_toolbox/gammatonegram")
	S_matlab = abs(specgram(s,nfft,fs,nwin,nwin-nhop));
	W_matlab = fft2gammatonemx(512, fs, 23, 1, 150, 4000, 512/2+1)
	(G_matlab, F_matlab) = gammatonegram(s, fs, 0.025, 0.010, 23, 150, 4000, 1)
end
@mget S_matlab
@mget G_matlab
@mget F_matlab
@mget W_matlab
#close_default_msession()

@test_approx_eq W W_matlab
@test_approx_eq F F_matlab'
@test_approx_eq G G_matlab
