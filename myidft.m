function [ f t ] = myidft( F,hz )
% Inverse Discrete Fourier Transform (not fast!)
% Bryan Kaiser
% 7/14/15

% This function computes the inverse Fourier 
% transform of the transformed signal "F" for 
% a sampling frequency in Hertz "hz" back into 
% a timeseries of uniform timestep "dt."

% f = the signal vector in the time-domain
% F = the signal vector in the frequency-domain
%        at the correct amplitude, including for the
%        dc mean value.
% dt = sample time step

% SUPPORTING FUNCTIONS: none

%============================

% Necessary parameters
N = length(F); % length of signal
Fs = hz(2).*N; % sampling rate
dt = (Fs.^(-1)); 
t = [0:(N-1)].*dt;

% Inverse Fourier transform (my way):
f = zeros(N,N);
F(1) = 2*F(1);
for k = 1:N % time loop
    for n = 1:N % frequency loop
        f(k,n) = F(n)*exp(sqrt(-1)*(2*pi*(n-1)/N)*(k-1)); 
    end
end
f = sum(f,2)./2;

end

