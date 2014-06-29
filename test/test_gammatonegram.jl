using Auditory, MATLAB, Base.Test, WAV

(s, fs) = wavread("/Users/jfsantos/Projects/SRMR_toolbox/gammatonegram/sa2.wav")
G, F = gammatonegram(s[:,1], fs, 0.025, 0.010, 23, 150, 4000, 1)

restart_default_msession()
@matlab begin
	addpath("/Users/jfsantos/Projects/SRMR_toolbox/gammatonegram")
    (s, fs) = wavread("/Users/jfsantos/Projects/SRMR_toolbox/gammatonegram/sa2.wav")
	(G_matlab, F_matlab) = gammatonegram(s, fs, 0.025, 0.010, 23, 150, 4000, 1)
end
@mget G_matlab
@mget F_matlab
close_default_msession()

@test_approx_eq G G_matlab
@test_approx_eq F F_matlab'
