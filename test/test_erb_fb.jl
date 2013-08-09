using auditory, MATLAB

# testing auditory filterbank design

restart_default_msession()

fb = make_erb_filterbank(16000, 23, 150)

@matlab begin
    addpath("/home/jfsantos/Dropbox/INRS/Summer 2013/Thesis/SRMR_normalization/auditory")
    fb_matlab = MakeERBFilters(16000, 23, 150)
end

@mget fb_matlab

# if max(abs(fb_matlab-fb)) < 1E-10
#     println("Implementation results under specified tolerance.")
# else
#     println("ERROR: Julia implementation not under specified tolerance.")
# end

# close_default_msession()

t = [1:32000]
x = sin(2*pi*440*t)
y = cochlear_filterbank(x, fb)

@mput x
@matlab y_matlab = ERBFilterBank(x, fb_matlab)
@mget y_matlab

sum(abs((y-y_matlab).^2))