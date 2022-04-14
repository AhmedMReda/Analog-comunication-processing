close all
clear
clc
Fs = 200;               % sampling frequency
Ts = 1/Fs;              % sampling time
resolution = 0.01;
T = (1/resolution)/2;
t = -T:Ts:T-Ts; % time vector
m_t = cos(2*pi*9*t);
%to make the cosine function to be defined form 0 to 3 and zero elsewhere
t03 = [zeros(1,10000) ones(1,600) zeros(1,9400)]; 
m_t = m_t.*t03;
n = length(m_t);
n2 = 2^nextpow2(n);
m_f = fft(m_t,n2)/Fs;
m_f0 = fftshift(m_f);
f = (-n2/2:n2/2-1)*(Fs/n2); % zero-centered frequency range
subplot(2,1,1);
plot(t,m_t);
axis([-1 4 -1 1]);
xlabel('Time (S)');
ylabel('Amplitude (V)');
title('Time Domain Signal');
subplot(2,1,2);
x=0:0.1:100; %the duaration of a line
g1=0.01*max(abs(m_f0)); %to make line at 1% of max to get the B.W 
plot(f,abs(m_f0),x,g1+x*0);
grid on
grid minor
axis([-20 20 0 1.5])
xlabel('Frequency (HZ)');
ylabel('Amplitude');
title('Frequency Domain Signal');
legend('Frequency Domain Signal','line at 1% of the max')
% % % % % % %analytical solution for m_f0% % % % % %
m_f0_A = 0.5*(((exp((1i*6*pi).*(9-f))-1)./((1i*2*pi).*(9-f)))+((exp((-1i*6*pi).*(9+f))-1)./((1i*2*pi).*(9+f))));
figure
plot(f,abs(m_f0),f,abs(m_f0_A));
axis([-20 20 0 1.5])
xlabel('Frequency (HZ)');
ylabel('Amplitude');
title('Frequency Domain Signal & the analytical solution');
legend('Frequency Domain Signal','the analytical solution')
% % % % % % %carrier% % % % % %
c_t = cos(2*pi*50*t);
c_f = fft(c_t,n2)/Fs;
c_f0 = fftshift(c_f);
% % % % % % %DSB-SC-point_a% % % % % %
s = m_t.*c_t;
s_f = fft(s,n2)/Fs;
s_f0 = fftshift(s_f);
figure
subplot(3,1,1);
g2=0.01*max(abs(s_f0)); %to make line at 1% of max to get the B.W 
plot(f,abs(s_f0),x,g2+x*0);
grid on
grid minor
axis([-100 100 0 0.8]);
xlabel('Frequency (HZ)');
ylabel('Amplitude');
title('DSB-SC modulated signal');
legend('DSB-SC modulated signal','line at 1% of the max')
% % % % % %DSB-SC-point_b% % % % % %
v = s.*c_t;
v = v.*2; %to remove the half
v_f = fft(v,n2)/Fs;
v_f0 = fftshift(v_f);
%design LPF (B.W = 20HZ) these numbers by ratios for ex:
% to make zeros from -20 to 20 (40) and graph form -100 to 100 (200) 
% and length(f)= 32768 then ones must be (40/200*32768 = 6554)
LPF = [zeros(1,13107) ones(1,6554) zeros(1,13107)]; %LPF (B.W = 20HZ)
m2_f0 = v_f0.*LPF;
m2_f = ifftshift(m2_f0);
m2_t = ifft(m2_f,n2)*Fs;
m2_t(20001:32768) = []; %remove the unnecessary data and to equalize the length
subplot(3,1,2);
plot(t,m2_t);
axis([-1 4 -1 1]);
xlabel('Time (S)');
ylabel('Amplitude (V)');
title('DSB-SC demodulated Signal');
subplot(3,1,3);
plot(t,m2_t,t,m_t);
legend('DSB-SC demodulated Signal','original signal')
axis([-1 4 -1 1]);
xlabel('Time (S)');
ylabel('Amplitude (V)');
title('DSB-SC demodulated Signal & original signal');
% % % % % %%SSB-point_a% % % % % %
%design BPF (B.W = 20HZ,center at 60HZ)
BPF = [zeros(1,4915) ones(1,3277) zeros(1,16384) ones(1,3277) zeros(1,4915)]; 
s2_f0 = s_f0.*BPF;
figure
subplot(3,1,1);
g3=0.01*max(abs(s2_f0)); %to make line at 1% of max to get the B.W 
plot(f,abs(s2_f0),x,g3+x*0);
grid on
grid minor
axis([-100 100 0 0.8]);
xlabel('Frequency (HZ)');
ylabel('Amplitude');
title('USB modulated signal');
legend('USB modulated signal','line at 1% of the max')
% % % % % %%SSB-point_b% % % % % %
v2_f0 = conv(s2_f0,c_f0,'same');
m3_f0 = v2_f0.*LPF;
m3_f = ifftshift(m3_f0);
m3_t = ifft(m3_f,n2)*Fs;
m3_t(20001:32768) = []; %remove the unnecessary data and to equalize the length
m3_t = (m3_t/max(-m3_t)); % correct the amp
subplot(3,1,2);
plot(t,m3_t);
axis([-1 4 -1 1])
xlabel('Time (S)');
ylabel('Amplitude (V)');
title('USB demodulated Signal');
subplot(3,1,3);
plot(t,m3_t,t,m_t);
legend('USB demodulated Signal','original signal')
axis([-1 4 -1 1])
xlabel('Time (S)');
ylabel('Amplitude (V)');
title('USB demodulated Signal & original signal');
figure
plot(f,LPF,f,BPF);
legend('LPF','BPF')
grid on
% % % % % %%AM% % % % % %
Fs=200; % sampling freq
Ts=1/(Fs);
% set no of samples
resolution = 0.01;
T = (1/resolution)/2;
t = -T:Ts:T-Ts; % time vector
mt = cos(2*pi*9*t);
%to make the cosine function to be defined form 0 to 3 and zero elsewhere
t03 = [zeros(1,10000) ones(1,600) zeros(1,9400)];
mt = mt.*t03; % basic signal
n = length(mt);
n2 = 2^nextpow2(n);
f = (-n2/2:n2/2-1)*(Fs/n2); % zero-centered frequency range
%parameters
Am=1;
Ac=2.5;
fm=9;
fc=50;
wc1=2*pi*fm;
wc2=2*pi*fc;
k=.4;
% carier signal
ct=Ac*cos(wc2*t);
%-------------------------------------------------------------
% modulated signal
s=(1+k*mt).*ct;
%-------------------------------------------------------------
% spectrum
S = fft(s, n2)/Fs;
S_shift=fftshift(abs(S));
subplot(2,1,1);
plot( t, s);
axis ([0 3 -5 5]);
title ('AM signal');
xlabel('Time');
ylabel('amplitude');
%----------------------------------------------------------------
subplot(2,1,2);
x=0:0.1:100; %the duaration of a line
g=0.015; %to make line at 1% of the signal max to get the B.W
plot(f,abs(S_shift),x,g+x*0);
axis ([-100 100 0 2]);
title('spetrum of AM');
xlabel('freq.');
ylabel('amplitude');
legend('spetrum of AM','line at 1% of the max')
grid on
grid minor