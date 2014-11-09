immutable ModulationFilterbank{T<:Filter}
	filters::Vector{T}
end

function filt(fb::ModulationFilterbank, x)
	output = zeros(length(x), length(fb.filters))
	@inbounds begin
		for k=1:length(fb.filters)
			output[:,k] = filt(fb.filters[k], x)
		end
	end
	return output
end

function make_modulation_filter(w0, Q)
    W0 = tan(w0/2)
    B0 = W0/Q
    b = [B0, 0, -B0]
    a = [(1 + B0 + W0^2), (2*W0^2 - 2), (1 - B0 + W0^2)]
    BiquadFilter(b[1], b[2], b[3], a[1], a[2], a[3])
end

function modulation_filterbank(mf, fs, q)
    ModulationFilterbank([make_modulation_filter(w0, q) for w0 in 2*pi*mf/fs])
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
