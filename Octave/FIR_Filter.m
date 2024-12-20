close all;
clear all;
clf;

f1 = 10000;
f2 = 15000;
delta_f = f2-f1;
Fs = 192000;
dB  = 40;
N = dB*Fs/(22*delta_f);

f =  [f1 ]/(Fs/2)
hc = fir1(round(N)-1, f,'low')

figure
plot((-0.5:1/4096:0.5-1/4096)*Fs,20*log10(abs(fftshift(fft(hc,4096)))))
axis([0 20000 -60 20])
title('Filter Frequency Response')
grid on