function [ F hz ] = mydft( f,dt )
% Discrete Fourier Transform (not fast!)
% Bryan Kaiser
% 7/14/15

% This function computes the discrete Fourier 
% transform of the signal "f" for a uniform time
% step "dt".

% Note: the output is not normalized, 
% and it is complex. Not: |F(omega)|

% f = the signal vector in the time-domain
% F = the signal vector in the frequency-domain
%        at the correct amplitude, including for the
%        dc mean value.
% dt = sample time step

% SUPPORTING FUNCTIONS: none

%============================

% Necessary parameters
N = length(f); % length of signal
%bnds = [1:floor(N/2)+1]; % single side of spectrum
hz = (0:N-1); % Hz
Fs = 1/dt; % 1/s

% Fourier transform into the frequency domain:
F = zeros(N,N);
for n = 1:N % frequency loop
    for k = 1:N % time loop
        F(n,k) = f(k)*exp(-sqrt(-1)*(2*pi*(n-1)/N)*(k-1)); 
    end
end
F = sum(F,2);
%F = (abs(F)).*2/N; % frequency domain signal 
% (single side), correct magnitude
F(1) = F(1)/2; % dc correction
hz = hz.*Fs/N; % Hz

end %---------------------------------------------------------

% To plot the output:

% % Compare to mydft.m function:
% [ F hz ] = mydft( f,dt );
% %bnds = [1:floor(N/2)+1]; % single side of spectrum
% figure
% if Fs >= 100
% semilogx(hz,F,'b'); else
% plot(hz,F,'b');  
% end
% axis([-1,(N/2+20)*Fs/N,0,8])
% xlabel('f (Hz)'); ylabel('|F[n]|')
% title('Conventional Spectrum')
