module auditory
export hilbert, erb_space, make_erb_filterbank, erb_filterbank, compute_modulation_cfs, make_modulation_filter, modulation_filterbank

function hilbert(x::Array{Float64})
# Return the Hilbert transform of x (a real signal).
# Code inspired by Scipy's implementation, which is under BSD license.
    X = vcat(rfft(x), zeros(int(floor(length(x)/2))-1))
    N = length(X)
    h = zeros(N)
    if N % 2 == 0
        h[1] = h[N/2+1] = 1
        h[2:N/2] = 2
    else
        h[1] = 1
        h[2:(N+1)/2+1] = 2
    end
    return ifft(X.*h)
end

function erb_filterbank(x, fcoefs)
    A0  = fcoefs[:,1]
    A11 = fcoefs[:,2]
    A12 = fcoefs[:,3]
    A13 = fcoefs[:,4]
    A14 = fcoefs[:,5]
    A2  = fcoefs[:,6]
    B0  = fcoefs[:,7]
    B1  = fcoefs[:,8]
    B2  = fcoefs[:,9]
    gain= fcoefs[:,10]	

    output = zeros(length(x), size(gain,1))
    
    for chan = 1:size(gain,1)

	y1=filt([A0[chan]/gain[chan], A11[chan]/gain[chan], A2[chan]/gain[chan]], [B0[chan], B1[chan], B2[chan]], x)
	y2=filt([A0[chan], A12[chan], A2[chan]], [B0[chan], B1[chan], B2[chan]], y1)
	y3=filt([A0[chan], A13[chan], A2[chan]], [B0[chan], B1[chan], B2[chan]], y2)
	y4=filt([A0[chan], A14[chan], A2[chan]], [B0[chan], B1[chan], B2[chan]], y3)
	output[:, chan] = y4
    end
    return output
end

function make_erb_filterbank(fs, num_channels, low_freq, EarQ = 9.26449, minBW = 24.7, order = 1)
    T = 1/fs
    if length(num_channels) == 1
	cf = erb_space(low_freq, fs/2, num_channels)
    else
	cf = num_channels
	if size(cf,2) > size(cf,1)
	    cf = cf'
	end
    end
    ERB = ((cf/EarQ).^order + minBW^order).^(1/order)
    B = 1.019*2*pi*ERB
    A0 = T
    A2 = 0
    B0 = 1
    B1 = -2*cos(2*cf*pi*T)./exp(B*T)
    B2 = exp(-2*B*T)

    A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2
    A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2
    A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2
    A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./exp(B*T))/2

    gain = abs((-2*exp(4*im*cf*pi*T)*T + 2*exp(-(B*T) +
      2*im*cf*pi*T).*T.*(cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))*
      sin(2*cf*pi*T))) .* (-2*exp(4*im*cf*pi*T)*T + 2*exp(-(B*T) +
      2*im*cf*pi*T).*T.* (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) *
      sin(2*cf*pi*T))).* (-2*exp(4*im*cf*pi*T)*T + 2*exp(-(B*T) +
      2*im*cf*pi*T).*T.* (cos(2*cf*pi*T) - sqrt(3 +
      2^(3/2))*sin(2*cf*pi*T))) .* (-2*exp(4*im*cf*pi*T)*T +
      2*exp(-(B*T) + 2*im*cf*pi*T).*T.* (cos(2*cf*pi*T) + sqrt(3 +
      2^(3/2))*sin(2*cf*pi*T))) ./ (-2 ./ exp(2*B*T) -
      2*exp(4*im*cf*pi*T) + 2*(1 + exp(4*im*cf*pi*T))./exp(B*T)).^4)
    
    allfilts = ones(length(cf),1);
    return [A0*allfilts A11 A12 A13 A14 A2*allfilts B0*allfilts B1 B2 gain], ERB;
end

function erb_space(low_freq, high_freq, num_channels, EarQ = 9.26449, minBW = 24.7, order = 1)
    # All of the following expressions are derived in Apple TR #35, "An
    # Efficient Implementation of the Patterson-Holdsworth Cochlear
    # Filter Bank."  See pages 33-34.
    cfArray = -(EarQ*minBW) + exp([1:num_channels]*(-log(high_freq + EarQ*minBW) + log(low_freq + EarQ*minBW))/num_channels) * (high_freq + EarQ*minBW)
end

function make_modulation_filter(w0, Q)
    W0 = tan(w0/2)
    B0 = W0/Q
    b = [B0, 0, -B0]
    a = [(1 + B0 + W0^2), (2*W0^2 - 2), (1 - B0 + W0^2)]
    b = b/a[1]
    a = a/a[1]
    return b, a
end

function modulation_filterbank(x, mf, fs, q)
    N = length(mf)
    out = zeros(length(x),N)
    for k=1:N
        w0 = 2*pi*mf[k]/fs
        (b3,a3) = make_modulation_filter(w0,q)
        out[:,k] = filt(b3, a3, x)
    end
    return out
end

function compute_modulation_cfs(min_cf, max_cf, n)
    spacing_factor = (max_cf/min_cf)^(1/(n-1))
    cfs = zeros(n)
    cfs[1] = min_cf
    for k=2:n
        cfs[k] = cfs[k-1]*spacing_factor
    end
    return cfs
end

end #module