% Jack, Keith. Video Demystified (p. 408). Elsevier Science. Kindle Edition. 

close all;
clear all;
clc;

pkg load signal

% Sample clock rate (13.5 Mhz)
Fs = 13.5e6;

% Line frequency (15.625 Khz)
% PAL - 64ms line frequency
Fh = 15.625e3;

% Calculate HCOUNT
H = 1/(Fh/Fs);

% The colour subcarrier frequency is based on an odd multiple of one-quarter the line rate, using the factor 1135/4
% The 4.43361875 Mhz of the color carrier is the result of 283.75 color clock cycles per line
% plus a 25 Hz offset to avoid interferences.
% Since the line frequency (number of lines per second) is 15625 Hz (625 lines x 50 Hz / 2)

% Determine P1 from 2024
Fsc = 4.43361875e6;
F = Fsc/Fs;
q1 = 2048;
p1 = (Fsc * q1) / Fs;
p1 = floor(p1);

%0.328125


% Determine P2 for PAL(B,D,G,H,I,N)
x = 1135/4 + 1/625;
p2 = 8192 * x - 4 * H * p1;
p2 = floor(p2);

% Determine P3
p3 = (Fsc * 625) / Fs;

% Constraints
% Since [p] is of finite word length, the DTO output frequency can be varied only in steps. 
% With a [p] word length of [w], the lowest [p] step is 0.5w and the lowest DTO frequency step is:
w = 11 % word length
w_low = 0.5 * w;
Fstep = Fs/(2^w);

% Note that the output frequency cannot be  greater than half the input frequency. 
% This  means that the output frequency FSC can only  be varied by the increment [p] and within the  range: 
Frange = 0 < Fsc < Fs/2;

t = 0:1/(4048*Fsc):1/(2*Fsc);
sin_rom = sin(2*pi*Fsc*t);

xn = zeros(1, 2048);
xn(1) = p1;
for i = 2:2048
  xn(i) = mod(xn(i-1) + p1, 2048);
endfor

figure;
plot(sin_rom);
%ylim([0 2048])
