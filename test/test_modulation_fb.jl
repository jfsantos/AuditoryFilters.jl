require("auditory.jl")
using auditory, MATLAB

x = [1, zeros(1023)]
mf = compute_modulation_cfs(4,128,8)
y = modulation_filterbank(x, mf, 16000, 2)

restart_default_msession()
@matlab begin
    addpath("/home/jfsantos/host/SRMRtoolbox/auditory");
    x = [1; zeros(1023,1)];
    mf_matlab = computeModulationCFs_mod(4, 128, 8);
    y_matlab = ModulationFilterbank(x, mf_matlab, 16000, 2);
end
@mget mf_matlab
@mget y_matlab

