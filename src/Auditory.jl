module Auditory

import DSP: Filter, TFFilter, BiquadFilter, SOSFilter, filt, spectrogram
import Base: convert
export erb_space, make_erb_filterbank, erb_filterbank, compute_modulation_cfs, make_modulation_filter, modulation_filterbank, gammatonegram

include("ERBFilterbank.jl")
include("ModulationFilterbank.jl")
include("Gammatonegram.jl")

end #module