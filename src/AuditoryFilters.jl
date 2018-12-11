__precompile__()
module AuditoryFilters

using DSP, Compat
if VERSION >= v"0.7.0-DEV.602"
    import ..DSP: filt, filt!
else
    import Base: filt, filt!
end
import Base: convert
export erb_space, make_erb_filterbank, compute_modulation_cfs, make_modulation_filter,
       modulation_filterbank, gammatonegram, fft2gammatonemx,
       ERBFilterbank, ModulationFilterbank, filt

abstract type Filterbank end

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
