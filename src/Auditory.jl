module Auditory

using DSP
import Base: convert, filt
export erb_space, make_erb_filterbank, erb_filterbank, compute_modulation_cfs, make_modulation_filter, modulation_filterbank, gammatonegram, 
	   ERBFilterbank, ModulationFilterbank

include("ERBFilterbank.jl")
include("ModulationFilterbank.jl")
include("Gammatonegram.jl")

end #module