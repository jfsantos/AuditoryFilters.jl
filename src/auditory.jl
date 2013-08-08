module auditory
export hilbert

function hilbert(x)
# Return the Hilbert transform of x.
# Code inspired by Scipy's implementation, which is under BSD license.
    X = fft(x)
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

end #module