using Auditory, MATLAB, Base.Test

x = [1, zeros(1023)]
mf = compute_modulation_cfs(4,128,8)
mfb = modulation_filterbank(mf, 16000, 2)
y = filt(mfb, x)

restart_default_msession()
@matlab begin
	addpath("/Users/jfsantos/Projects/SRMR_toolbox/auditory")
    x = [1; zeros(1023,1)];
    mf_matlab = computeModulationCFs_mod(4, 128, 8);
    y_matlab = modulationFilterBank(x, mf_matlab, 16000, 2);
end
@mget mf_matlab
@mget y_matlab
close_default_msession()

@test_approx_eq y y_matlab'
@test_approx_eq mf mf_matlab
