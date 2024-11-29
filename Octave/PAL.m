close all;
clear all;
clc;

pkg load signal

eq_pulses = 5;
h_lines = 620;

% Generate horizontal sync
fs = 40e6;

% Horizontal sync pulses (4.7us)
e1 = 1280e-6;
t1 = 0:1/fs:e1;
f1 = 1/64e-6;
h = 0.15 * -square(2*pi*f1*t1, 7.34375) + 0.15;

% Equalizing pulses
e2 = 160e-6;
t2 = 0:1/fs:e2;
f2 = 1/32e-6;
ep1 = 0.15 * -square(2*pi*f2*t2, 7.34375) + 0.15;

% Equalizing pulses
e3 = 160e-6;
t3 = 0:1/fs:e3;
f3 = 1/32e-6;
ep2 = 0.15 * -square(2*pi*f3*t3, 85.3125) + 0.15;
ep = [ep2 ep1 h];

% Color burst
f4 = 4.43361875e6;
cb_start = 1285.6e-6/(1/fs);
cb_stop  = cb_start + (1/f4*10)/(1/fs);
t = 0;
for i = 1:columns(ep)
  % Sub-carrier
  if (i >= cb_start && i <= cb_stop)
    cs(i) = ep(i) + 0.15 * sin(2*pi*f4*t);
  else
    cs(i) = ep(i);
  endif
  t += 1/fs;
endfor

figure;
plot(cs); % o-r
ylim([0 1])