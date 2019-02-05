function [ CWT scale ] = morletCWT( f,t,F,omega )
% Morlet Continuous Wavelet Transform
% Bryan Kaiser
% 7/23/15

% This function takes a discrete signal "f" and its 
% discrete time stamps and generates the 
% continuous wavelet transform of the function 
% "f" over t. NOW IN PARALLEL!!

% SUPPORTING FUNCTIONS: mydft.m

%============================

dt = diff(t(1:2));
N = length(t);

% CWT Morlet wavelet scales ============
dj = 0.1; % scale spacing 
omega0 = 6; % nondimensional frequency chosen for admissibility condition
s0 = 2*dt;
J = floor(1/dj*log2((N*dt)/s0)); % Maximum possible number of scales
j = [0:ceil(J)]; % scale s indices
scale = zeros(1,length(j));
for s = 1:length(j)
    scale(s) = s0*2^(j(s)*dj); % seconds, continous scales 
end


% FFT of time series =================

% -- My way:
% [ F hz ] = mydft( f,dt ); % The output must be complex! F(\omega), not |F(omega)|
% omega = hz.*(2*pi); % Angular frequency
% for k = 1:N
%     if k > N/2
%         omega(k) = -omega(k);
%     end
% end

% % -- Matlab way:
% NFFT = 2^nextpow2(N); % Next power of 2 from length of y
% F = fft(f,NFFT)/N; % Fourier transform to frequency domain
% Fs = 1/dt;
% omega = Fs/2*linspace(0,1,NFFT/2+1); % Frequency domain
% for k = 1:N
%     if k > N/2
%         omega(k) = -omega(k);
%     end
% end

% Continuous Wavelet Transform =========
CWT = nan(length(t),length(scale));

% omega = zeros(1,N); % omega_k
% for k = 1:N
%     omega(k) = (2*pi*(k-1))/(N*dt);
%     %if k > N/2
%     %    omega(k) = -omega(k);
%     %end
% end

parpool('local',2);
for s = 1:length(scale) % for each scale (maximum N)   
    
    Ws0 = nan(1,N); % Morlet mother wavelet
    
    % Eqn (6)
    parfor k = 1:N
        if omega(k) > 0 % Heaviside step function
        Ws0(k) = ((2*pi*scale(s)/dt)^(0.5))*(pi^(-0.25))*exp(-((scale(s)*omega(k)-omega0)^2)/2); 
        else 
            Ws0(k) = 0;
        end % Morlet mother wavelet for scale s
    end 
     
   Ws = nan(N,N);

   % Eqn (4)
   parfor n = 1:N % time loop 
       for k = 1:N % frequency loop
        Ws(n,k) = F(k)*(conj((Ws0(k))))*exp(sqrt(-1)*omega(k)*t(n)); 
       end  % Complex conjugate of normalized daugther wavelet times the FT of the time 
   end % series, for scale s, prepped for IFT (note ts = [0,1,...,N-1]*dt)
   CWT(:,s) = (sum(Ws,2)); 
   % The WT for scale s, a time series.
   
   % Why does the complex conjugate not matter, same thing with
   % omega_k????????
   
 
end % for scale s  
delete(gcp);

%figure
%plot(1:N,omega)
% figure
% plot(omega,Ws0)
% Ws0(1)
% Ws0(end)
end 

