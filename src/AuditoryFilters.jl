module AuditoryFilters

using DSP
import Base: convert, filt
export erb_space, make_erb_filterbank, erb_filterbank, compute_modulation_cfs, make_modulation_filter, modulation_filterbank, gammatonegram, 
	   ERBFilterbank, ModulationFilterbank

abstract Filterbank

nchannels(f::Filterbank) = length(f.filters)

function filt(fb::Filterbank, x)
	output = zeros(length(x), length(fb.filters))
	@inbounds begin 
		for k = 1:length(fb.filters)
			output[:, k] = filt(fb.filters[k], x)
		end
	end
	return output
end

include("ERBFilterbank.jl")
include("ModulationFilterbank.jl")
include("Gammatonegram.jl")

end #module