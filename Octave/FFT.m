%FFT example
close all;
clear all;
clc;

Fs = 10e6;
Tdur = 100e-6;
FStart = 500e3;
FStop = 1500e3;

tvec = 0:1/Fs:Tdur-1/Fs;
xreal = sin(2*pi*(FStart*tvec + (FStop-FStart)/(2*Tdur)*tvec.^2));

figure
plot(tvec,xreal)
xlabel('Time in [sec]')
ylabel('Amplitude in [V]')
grid on

Nfft = length(xreal);
fvec = (-Nfft/2:Nfft/2-1)*Fs/Nfft;
Xreal = fftshift(fft(xreal,Nfft));

figure
plot(fvec,abs(Xreal))
xlabel('Frequency in [Hz]')
ylabel('Amplitude')
grid on