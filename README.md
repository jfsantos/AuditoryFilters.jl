[![Build Status](https://travis-ci.org/jfsantos/AuditoryFilters.jl.svg?branch=master)](https://travis-ci.org/jfsantos/AuditoryFilters.jl)
[![Coverage Status](https://coveralls.io/repos/jfsantos/AuditoryFilters.jl/badge.png)](https://coveralls.io/r/jfsantos/AuditoryFilters.jl)


This package provides auditory filter implementations in Julia. The following implementations are available:

- A gammatone filterbank based on Malcolm Slaney's [Auditory Toolbox](https://engineering.purdue.edu/~malcolm/interval/1998-010/);
- a gammatone-like spectrogram implementation, based on Dan Ellis' [implementation](http://www.ee.columbia.edu/ln/rosa/matlab/gammatonegram/);
- a modulation filterbank based on the paper titled "Characterizing frequency selectivity for envelope
fluctuations" (2000), by Ewert and Dau. 
