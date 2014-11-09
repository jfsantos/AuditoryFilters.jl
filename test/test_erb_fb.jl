# testing ERB filterbank design
fb = make_erb_filterbank(16000, 23, 150)
fb_matlab = readcsv(open(joinpath(dirname(@__FILE__), "data", "ERB_filter_coeffs.csv")))
@test length(fb.filters) == size(fb_matlab, 1)
for k=1:length(fb.filters)
	sos = fb.filters[k]
	biquads = sos.biquads
	g = sos.g
	for b=1:4
		@test_approx_eq biquads[b].b0 fb_matlab[k,1]
		@test_approx_eq biquads[b].b2 fb_matlab[k,6]
		@test_approx_eq biquads[b].a1 fb_matlab[k,8]
		@test_approx_eq biquads[b].a2 fb_matlab[k,9]
	end	 
	@test_approx_eq biquads[1].b1 fb_matlab[k,2]
	@test_approx_eq biquads[2].b1 fb_matlab[k,3]
	@test_approx_eq biquads[3].b1 fb_matlab[k,4]
	@test_approx_eq biquads[4].b1 fb_matlab[k,5]
	@test_approx_eq 1/g fb_matlab[k,10]
end

# testing ERB filterbank response
x = zeros(10000)
x[1] = 1.0
y = filt(fb, x)
y_matlab = readcsv(open(joinpath(dirname(@__FILE__), "data", "ERB_filter_response.csv")))
@test_approx_eq y y_matlab'
