%%% Performs forward Fourier Transform of function A, with user-set 
%%% resolution, R = (size(A)+2*padSize)

function [a] = FT(A,padSize)
    A = padarray(A,[padSize,padSize]);
    a = fftshift(fft2(A));  
    
end