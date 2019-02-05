% Morlet Wavelet Transform
% Bryan Kaiser
% 7/24/15

% This script generates a test signal and then 
% proceeds to plot the wavelet transform of the 
% signal, using a Morlet wavelet.

% INPUT: none

% OUTPUT: none

% SUPPORTING FUNCTIONS: morletCWT.m or 
% mycwt.m

% PROBLEM: Why does the CWT signal only 
% extend to t=300 even though in reality it 
% extends to 500? 
% ANSWER: Probably because a single-sided 
% fft.m was used!!

%============================

close all; clear all; clc;

% % Signal ========================
% N = 1E3; % total number of samples (multiple of 2)
% T = 1E2;
% Fs = N/T; % sampling rate
% dt = 1/Fs; % time step
% ts = zeros(1,N); % time samples
% x = zeros(1,N); % sampled signal
% for i = 1:N
%     ts(i) = dt*(i-1); % t = [0,dt,2*dt,...,(N-1)*dt]
%     if i < N/2
%         x(i) = sin(0.1*pi*ts(i)); % 1/10 Hz
%     else
%         x(i) = sin(pi*ts(i)); % 1 Hz
%     end    
% end
% T = ts(end); % final time
% % THIS SIGNAL IS IDENTICAL TO THE PYTHON ONE!

% % Signal ========================
% N = 2E3; % total number of samples (multiple of 2)
% T = 2E2;
% Fs = N/T; % sampling rate
% dt = 1/Fs; % time step
% ts = zeros(1,N); % time samples
% x = zeros(1,N); % sampled signal
% for i = 1:N
%     ts(i) = dt*(i-1); % t = [0,dt,2*dt,...,(N-1)*dt]
%     if i < N/2
%         x(i) = sin(0.1*pi*ts(i)); % 1/10 Hz
%     else
%         x(i) = sin(pi*ts(i)); % 1 Hz
%     end    
% end
% T = ts(end); % final time
% % wall time 7min for N = 2000

% % Signal ========================
% N = 5000; % total number of samples (multiple of 2)
% T = 500;
% Fs = N/T; % sampling rate
% dt = 1/Fs; % time step
% ts = zeros(1,N); % time samples
% x = zeros(1,N); % sampled signal
% for i = 1:N
%     ts(i) = dt*(i-1); % t = [0,dt,2*dt,...,(N-1)*dt]
%     if i < N/2
%         x(i) = sin(0.1*pi*ts(i)); % 1/10 Hz
%     else
%         x(i) = sin(pi*ts(i)); % 1 Hz
%     end    
% end
% T = ts(end); % final time
% % wall time 1.1667hrs for N = 5000

% % Signal ========================
% N = 3E3; % total number of samples (multiple of 2)
% T = 2E3;
% Fs = N/T; % sampling rate
% dt = 1/Fs; % time step
% ts = zeros(1,N); % time samples
% x = zeros(1,N); % sampled signal
% for i = 1:N
%     ts(i) = dt*(i-1); % t = [0,dt,2*dt,...,(N-1)*dt]
%     if i < 2*N/3 & i >= N/3
%         x(i) = sin(0.005*2*pi*ts(i)); % 1/200 Hz, tau = 200s (CORRECT Hz!)
%     elseif i < N/3
%         x(i) = sin(0.025*2*pi*ts(i)); % 1/40 Hz, tau = 40s (CORRECT Hz!)
%     else
%         x(i) = sin(0.0025*2*pi*ts(i)); % 1/400 Hz, tau = 400s (CORRECT Hz!)
%     end    
% end
% T = ts(end); % final time
% % % wall time 18min for N = 3000
% % cwt_072415

% Signal ========================
N = 5E3; % total number of samples (multiple of 2)
T = 2500;
Fs = N/T; % sampling rate
dt = 1/Fs; % time step
ts = zeros(1,N); % time samples
x = zeros(1,N); % sampled signal
for i = 1:N
    ts(i) = dt*(i-1); % t = [0,dt,2*dt,...,(N-1)*dt]
    if i < 3*N/5 & i >= 2*N/5
        x(i) = sin(0.005*2*pi*ts(i)); % 1/200 Hz, tau = 200s (CORRECT Hz!)
    elseif i < 2*N/5
        x(i) = sin(0.025*2*pi*ts(i)); % 1/40 Hz, tau = 40s (CORRECT Hz!)
    else
        x(i) = sin(0.0025*2*pi*ts(i)); % 1/400 Hz, tau = 400s (CORRECT Hz!)
    end    
end
T = ts(end); % final time
% wall time 65min for N = 5000
 % cwt_072415b

% % Signal ========================
% load freqbrk; x = freqbrk; %figure; plot(x,'r');
% dt = 1; ts = [0:1:length(x)-1]; %hdf5write('signal.h5','/f',signal);
% T = ts(end);


% Signal Plot =====================
figure; plot(ts,x,'b','LineWidth',1.1); set(gca,'FontSize',12);
xlabel('$${t}\hspace{1mm}[s]$$','interpreter','latex','FontSize',14);
ylabel('$$f(t)$$','interpreter','latex','FontSize',14);
%ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') - [0.5 0 0]);
title('$$Signal$$','interpreter','latex','FontSize',14);

% Fourier Transform of signal ===========
NFFT = 2^nextpow2(N); % Next power of 2 from length of y
X = fft(x,NFFT)/N; % Fourier transform to frequency domain
hz = Fs*linspace(0,1,NFFT/2+1); % Frequency domain
% Plot single-sided amplitude spectrum.
figure; semilogx(hz,2*abs(X(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of x(t), fft.m')
xlabel('Frequency (Hz)'); ylabel('|X|')
%axis([-1,(N/2+20)*Fs/N,0,0.6])



% Continuous Wavelet Transform (Morlet) ====

%[ WT scale ] = mycwt( x,ts );

[ X hz ] = mydft( x,dt );
omega = hz.*(2*pi); % Angular frequency
for k = 1:N
    if k > N/2
        omega(k) = -omega(k);
    end
end
tic
[ WT scale ] = morletCWT( x,ts,X,omega );
CWTtime = toc 

% Wavelet Transform Variance ===========
WTr = real(WT);
WTi = imag(WT);
WTp = abs(WT).^2;

[N J] = size(WT);
dj = 0.1; % My sample rate was 10 hz (period 0.1)
% Could that cause problems?

Cd = 0.776 % Test: Morlet
% Cdj = zeros(N,J);
% for j = 1:J
%     Cdj(j) = WTr(1,j)/(scale(j)^(1/2));
% end
% Cd =  (dj*sqrt(dt))/(pi^(1/4))*sum(Cdj)

sigma2 = zeros(N,J);
for n = 1:N
    for j = 1:J
        sigma2(n,j) = WTp(n,j)/scale(j);
    end
end
sigma2 = (dj*sqrt(dt))/(Cd*N)*sum(sum(sigma2));

% Cone of Influence =================

% ADD!!!!!!!!!

% Plot of CWT =====================

N = length(x);
[tgrid,sgrid]=meshgrid(ts,scale);

% Real component
figure; set(gca,'FontSize',12); h2 = surf(tgrid,sgrid,WTr'); 
minima = min(min(WTr)); maxima = max(max(WTr)); 
colorbar; caxis([minima,maxima]); view(2); shading flat; grid off; 
axis tight; xlabel('t'); ylabel('scale')

% Complex component
figure; set(gca,'FontSize',12); h2 = surf(tgrid,sgrid,WTi'); 
minima = min(min(WTi)); maxima = max(max(WTi)); 
colorbar; caxis([minima,maxima]); view(2); shading flat; grid off; 
axis tight; xlabel('t'); ylabel('scale')

% Wavelet Power Spectrum (unnormalized)
figure; set(gca,'FontSize',12); h2 = surf(tgrid,sgrid,WTp'./sigma2); 
colorbar;  view(2); shading flat; grid off; 
axis tight; %caxis([10,100]); axis([0,T,scale(1),400]);
%colormap winter;
xlabel('$${t}\hspace{2mm}[s]$$','interpreter','latex','FontSize',12);
ylabel('$$scale\hspace{2mm}[s]$$','interpreter','latex','FontSize',12);
%ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') - [0.7 0 0]);
title('$$Wavelet\hspace{1mm}Power\hspace{1mm}Spectrum\hspace{1mm}|W_n(s)|^2/\sigma^2$$','interpreter','latex','FontSize',14);
%axis 'ij'; axis([0,205,0,4500]);       
%hold on; plot(ts,COI,'k','LineWidth',2);


