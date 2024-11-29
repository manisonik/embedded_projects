close all;
clear all;
clc;

pkg load signal

n = 5;
f = 6e6;

[zb,pb,kb] = butter(n,2*pi*f,'s');
[bb,ab] = zp2tf(zb,pb,kb);
w = logspace(1,12,5120); % Frequency vector
h = freqs(bb,ab,w);
mag = 20*log10(abs(h)); %// Magnitude in dB
pha = (180/pi)*angle(h); %// Phase in degrees

figure;
subplot(2,1,1);
semilogx(w, mag, 'r.');
ylim([-20 5])
xlabel('Frequency');
ylabel('Magnitude (dB)');
grid;

