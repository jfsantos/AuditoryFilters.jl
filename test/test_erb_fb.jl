using Auditory, MATLAB, Base.Test

# testing auditory filterbank design

restart_default_msession()

fb = make_erb_filterbank(16000, 23, 150)

@matlab begin
    addpath("/Users/jfsantos/Projects/SRMR_toolbox/auditory")
    fb_matlab = MakeERBFilters(16000, 23, 150)
end

@mget fb_matlab

x = zeros(1000)
x[1] = 1.0
y = filt(fb, x)

@mput x
@matlab y_matlab = ERBFilterBank(x, fb_matlab)
@mget y_matlab

close_default_msession()

@test_approx_eq y y_matlab'