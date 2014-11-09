# Converts an FFT into a gammatonegram. Based on the gammatonegram/fft2gammatonemx MATLAB 
# functions by Dan Ellis.
# D. P. W. Ellis (2009). "Gammatone-like spectrograms", web resource.
# http://www.ee.columbia.edu/~dpwe/resources/matlab/gammatonegram/

type Gammatonegram{T, F<:Real} <: DSP.Periodograms.TFR{T}
	amplitude::Matrix{T}
	frequencies::Vector{F}
	time::FloatRange{Float64}
end

function gammatonegram(x,sr::Integer,twin::Real,thop::Real,N::Integer,fmin,fmax,width)
	nfft = int(2^(ceil(log(2*twin*sr)/log(2))))
    nhop = iround(thop*sr)
    nwin = iround(twin*sr)
    (W,F) = fft2gammatonemx(nfft, sr, N, width, fmin, fmax)
    # perform FFT and weighting in amplitude domain
	S = stft(x, nwin, nwin-nhop; nfft=nfft, fs=sr, window=hanning)
	Y = 1/nfft*W*abs(S)
	Gammatonegram(Y, vec(F), ((0:size(Y,2)-1)*(nhop)+nwin/2)/sr)
end

function fft2gammatonemx(nfft::Integer, sr::Integer, N::Integer, width, fmin, fmax)
	# Allocating matrix to store the filterbank weights
	W = zeros(N, int(nfft/2+1));
	# Gammatone filterbank constants
	EarQ = 9.26449
	minBW = 24.7
	order = 1

	cfreqs = -(EarQ*minBW) + exp((1:N)'*(-log(fmax + EarQ*minBW) + log(fmin + EarQ*minBW))/N) * (fmax + EarQ*minBW);
	cfreqs = fliplr(cfreqs)
	GTord = 4
	ucirc = exp(im*2*pi*[0:(nfft/2)]/nfft)
	for k=1:N
	    cf = cfreqs[k];
	    ERB = width*((cf/EarQ).^order + minBW^order).^(1/order);
	    B = 1.019*2*pi*ERB;
	    r = exp(-B/sr);
	    theta = 2*pi*cf/sr;
	    pole = r*exp(im*theta);
	    T = 1/sr;
	    
		A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*
	                                                      cf*pi*T)./exp(B*T))/2; 
	    A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*
	                                                      cf*pi*T)./exp(B*T))/2;
	    A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*
	                                                      cf*pi*T)./exp(B*T))/2; 
	    A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*
	                                                      cf*pi*T)./exp(B*T))/2; 
	    zros = -[A11 A12 A13 A14]/T;
    
	    gain =  abs((-2*exp(4*im*cf*pi*T)*T + 
	                2*exp(-(B*T) + 2*im*cf*pi*T).*T.* 
	                (cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* 
	                 sin(2*cf*pi*T))) .* 
	               (-2*exp(4*im*cf*pi*T)*T + 
	                2*exp(-(B*T) + 2*im*cf*pi*T).*T.* 
	                (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * 
	                 sin(2*cf*pi*T))).* 
	               (-2*exp(4*im*cf*pi*T)*T + 
	                2*exp(-(B*T) + 2*im*cf*pi*T).*T.* 
	                (cos(2*cf*pi*T) - 
	                 sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* 
	               (-2*exp(4*im*cf*pi*T)*T + 2*exp(-(B*T) + 2*im*cf*pi*T).*T.* 
	                (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ 
	               (-2 ./ exp(2*B*T) - 2*exp(4*im*cf*pi*T) +  
	                2*(1 + exp(4*im*cf*pi*T))./exp(B*T)).^4);
	    W[k,:] = ((T^4)/gain) *
	        abs(ucirc-zros[1]).*abs(ucirc-zros[2]) .*
	        abs(ucirc-zros[3]).*abs(ucirc-zros[4]) .*
	        (abs((pole-ucirc).*(pole'-ucirc)).^-GTord);
	end
	W, cfreqs	
end
