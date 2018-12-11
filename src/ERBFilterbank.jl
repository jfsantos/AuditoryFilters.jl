struct ERBFilterbank{C,G,T<:Real,U<:Real,V<:Real} <: Filterbank
    filters::Vector{SecondOrderSections{C,G}}
    ERB::Vector{T}
    center_frequencies::Vector{U}
    fs::V
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
    ERB = ((cf/EarQ).^order .+ minBW.^order).^(1/order)
    B = 1.019*2*pi*ERB
    B0 = T
    B2 = 0.0
    A0 = 1.0
    A1 = -2*cos.(2*cf*pi*T)./exp.(B*T)
    A2 = exp.(-2*B*T)

    B11 = -(2*T*cos.(2*cf*pi*T)./exp.(B*T) .+ 2*sqrt(3+2^1.5)*T*sin.(2*cf*pi*T)./exp.(B*T))/2
    B12 = -(2*T*cos.(2*cf*pi*T)./exp.(B*T) .- 2*sqrt(3+2^1.5)*T*sin.(2*cf*pi*T)./exp.(B*T))/2
    B13 = -(2*T*cos.(2*cf*pi*T)./exp.(B*T) .+ 2*sqrt(3-2^1.5)*T*sin.(2*cf*pi*T)./exp.(B*T))/2
    B14 = -(2*T*cos.(2*cf*pi*T)./exp.(B*T) .- 2*sqrt(3-2^1.5)*T*sin.(2*cf*pi*T)./exp.(B*T))/2

    gain = abs.((-2*exp.(4*im*cf*pi*T)*T .+ 2*exp.(-(B*T) .+
      2*im*cf*pi*T).*T.*(cos.(2*cf*pi*T) .- sqrt(3 - 2^(3/2))*
      sin.(2*cf*pi*T))) .* (-2*exp.(4*im*cf*pi*T)*T .+ 2*exp.(-(B*T) .+
      2*im*cf*pi*T).*T.* (cos.(2*cf*pi*T) .+ sqrt(3 - 2^(3/2)) *
      sin.(2*cf*pi*T))).* (-2*exp.(4*im*cf*pi*T)*T .+ 2*exp.(-(B*T) .+
      2*im*cf*pi*T).*T.* (cos.(2*cf*pi*T) - sqrt(3 +
      2^(3/2))*sin.(2*cf*pi*T))) .* (-2*exp.(4*im*cf*pi*T)*T .+
      2*exp.(-(B*T) + 2*im*cf*pi*T).*T.* (cos.(2*cf*pi*T) + sqrt(3 +
      2^(3/2))*sin.(2*cf*pi*T))) ./ (-2 ./ exp.(2*B*T) -
      2*exp.(4*im*cf*pi*T) + 2*(1 .+ exp.(4*im*cf*pi*T))./exp.(B*T)).^4)

    C = typeof(B0)
    filters = Array{SOSFilter{C,C}}(undef,num_channels)
    for ch=1:num_channels
        biquads = Array{BiquadFilter{C}}(undef,4)
        biquads[1] = BiquadFilter(B0, B11[ch], B2, A0, A1[ch], A2[ch])
        biquads[2] = BiquadFilter(B0, B12[ch], B2, A0, A1[ch], A2[ch])
        biquads[3] = BiquadFilter(B0, B13[ch], B2, A0, A1[ch], A2[ch])
        biquads[4] = BiquadFilter(B0, B14[ch], B2, A0, A1[ch], A2[ch])
        filters[ch] = SOSFilter(biquads, 1/gain[ch])
    end
    ERBFilterbank(filters, ERB, cf, fs)
end

function erb_space(low_freq, high_freq, num_channels, EarQ = 9.26449, minBW = 24.7, order = 1)
    # All of the following expressions are derived in Apple TR #35, "An
    # Efficient Implementation of the Patterson-Holdsworth Cochlear
    # Filter Bank."  See pages 33-34.
    space = (1:num_channels).*(-log(high_freq + EarQ*minBW) + log(low_freq + EarQ*minBW))/num_channels
    cfArray = -(EarQ*minBW) .+ exp.(space) * (high_freq + EarQ*minBW)
end
