%function filtered_data = BPF(signal)



%%--- Program to design bandpass filter with basic mathematical equations
%%--- And will be helpful for those who dont have signal processing toolbox
%%% Created by:Ashwini Deshpande, (vd.ashwini@yahoo.co.in) 10/28/09

% Intial settings
Fs=.5;        % Sampling Frequency
t=(0:195)/Fs';   % Time Vector
L=196;         % FFT length
NFFT=256;      
f=Fs*linspace(0,1,NFFT); % Frequency Vector

%-- Generate composite signal with 100Hz and 200Hz frequency components --%
%comp_sig = signal;
comp_sig=sin(2*pi*.05*t)+sin(2*pi*.7*t);

%-- Calculate FFT of composite signal; --%
comp_sig_fft=fft(comp_sig,NFFT)/L;

bw=.045*pi;      %Bandwisth
fc=.055*pi;  %Center Frequency
%n=1:L;
%Compute Hamming window
for n=1:L
    hamm(n)=(0.54-0.46*cos(2*pi*n/L));
end
%Compute Filter
hsuup=(-(L-1)/2:(L-1)/2);
hideal1=hamm.*(2*(fc+bw)*(sin(2*(fc+bw)*hsuup/Fs)./(2*(2*fc+bw)*hsuup/Fs)));
hideal2=hamm.*(2*(fc-bw)*(sin(2*(fc-bw)*hsuup/Fs)./(2*(2*fc+bw)*hsuup/Fs)));
h_bpf=(hideal1-hideal2);
% Filtering in freq domain
h_bpf_fft=fft(h_bpf,NFFT)/L;
sfiltered_fft=comp_sig_fft.*h_bpf_fft;
filtered_data=real(ifft(sfiltered_fft));



%-- plot all signals ---%
figure;
subplot(411);
plot(t,comp_sig);grid on
%axis([0 0.1 -2 2])
title('Composite signal')

subplot(412)
plot(f,2*abs(comp_sig_fft));grid on;
%axis([0 5000 0 1.5])
title('Frequency components of composite signal')

subplot(413);
plot(t,filtered_data(1:L));hold on
title('After Filtering')

subplot(414);
plot(f,2*abs(sfiltered_fft));hold on;
%axis([0 500 0 10])
title('After filtering in frequency domain')


