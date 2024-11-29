close all;
clear all;
clc;

pkg load signal

% Start from nothing
clear;

% The sampling frequency in Hz.
fileID = fopen("SDS00004.csv");
[x] = dlmread (fileID, ";", 1, 0);
fclose(fileID);

% Clamp the signal
input = x(1:end,2);
min = min(input);
input = input - min;

% The sampling frequency in Hz (10 ns)
Fsam  = 1e9;

% Nyquist frequency in Hz.
% The Nyquist frequency is half your sampling frequency.
Fnyq = Fsam/2;

% The cut-off frequency of your Low pass filter in Hz.
% This frequency must be greater than 0 and less than Fnyq.
% Set the center frequency to 6 Mhz
Fc=6e6;

% Step 1: Lowpass filter at 6Mhz to remove any high-frequency components
% that might result in aliasing.
% Create a first-order Butterworth low pass at 6Mhz to filter out 
% The returned vectors are of legth n.
% Thus a first order filter is created with n = 2.
[b,a]=butter(2, Fc/Fnyq);
y=filtfilt(b,a,input);

% Apply the filter to the input signal and plot input and output.
z = fir1(500, 5e5/Fnyq, 'low');
p=filtfilt(z,1,y);

% Detect rise/fall of square wave
for i =  1:rows(p)
  if p(i) < 0.1
    z(i) = 0;
  else
    z(i) = 1;
  endif
endfor
z = diff(z);

% Sample the front porch for 3 us
st = 3e-6/1e-9;
d = 1e-6/1e-9;
avg = 0;
count = 0;
for i = 1:columns(z)
  if z(1,i) > 0.5
    c = i+d;
    s = p(c:c+st);
    avg = avg + mean(s);
    count = count + 1;
  endif
endfor
avg = avg / count;

% The front porch reference voltage of a PAL signal is 0.3V
% This is to simulate the AD8337
gain = 0.3/avg;
y = gain * y;
vgain = (gain - 12.65) / 19.7;

figure;
subplot (4, 1, 1)
plot(y)
subplot (4, 1, 2)
plot(p)
subplot (4, 1, 3)
plot(z)
subplot (4, 1, 4)
plot(s)
%freqz(b,a,[],Fsam);

% Remove 4.43361875e6 signal from PAL using a band-reject filter
%[b1,a1]=butter(2, 3e6/Fnyq);
%[b2,a2]=butter(2, 3e8/Fnyq, "high");
%b = conv(b1, b2);
%a = conv(a1, a2);
%freqz(b2,a2,4096,Fsam);

% Apply the filter to the input signal and plot input and output.
%Fc = 5e5;
%[b,a]=butter(1, Fc/Fnyq);

% Determine the blanking level by first lowpass filtering to about 0.5Mhz
%FcN = 4.43361875e6; % notch frequency
%fR = FcN/Fnyq;      % ratio of notch freq. to Nyquist freq.
%nW = 0.1;           % width of the notch

% Compute zeros
%nZ = [exp( sqrt(-1)*pi*fR ), exp( -sqrt(-1)*pi*fR )];

% Compute poles
%nP = (1-nW) * nZ;
%b = poly( nZ );   %  Get moving average filter coefficients
%a = poly( nP );   %  Get autoregressive filter coefficients


%% Design Notch filter for first frequency
%H(z) = (1 - 2 * cos(w1) * z^(-1) + z^(-2)) / (1 - 2 * p * cos(w1) * z^(-1) + p^2 * z^(-2))
%fs = 1e9;           % Sampling frequency [Hz]
%f1 = 4.43361875e6;  % Notch filter frequency 1 [Hz]
%p = 0.99; % The closer to one, the narrower the notch bandwidth
%w1 = 2*pi*f1 / fs; % Convert f1 to normalized frequency [+-2*pi]
%b = [1 , -2 * cos(w1), 1];      %IIR-Filter coefficients for numerator 1
%a = [1. -2 * p * cos(w1), p^2]; %IIR-Filter coefficients for denumerator 1

%% Create Filter object and plot
%freqz(b,a,[],fs);

% Create a first-order Butterworth low pass.
% The returned vectors are of legth n.
% Thus a first order filter is created with n = 2.
%[zb,pb,kb]=butter(1,2*pi*Fc,'s');
%[bb,ab] = zp2tf(zb,pb,kb);

%w = logspace(0,10); % Frequency vector
%h = freqs(bb,ab,w); %// Output of freqs
%mag = 20*log10(abs(h)); %// Magnitude in dB
%pha = (180/pi)*angle(h); %// Phase in degrees

%// Declare points
%wpt = [600, 7500];
%mpt = [20, -71];

%// Plot the magnitude as well as markers
%figure;
%subplot(2,1,1);
%semilogx(w, mag, wpt, mpt, 'r.');
%xlabel('Frequency');
%ylabel('Magnitude (dB)');
%grid;

%// Plot phase
%subplot(2,1,2);
%semilogx(w, pha);
%xlabel('Frequency');
%ylabel('Phase (Degrees)');
%grid;

%freqz(b,a, [1e6, 10e6],Fsam)
%subplot(2,1,1)
%ylim([-10 10])

% Apply the filter to the input signal and plot input and output.
%figure(2)
%output=filter(b,a,input);
%plot(output)

% Compute the fequency response
%[H,w] = freqz(bb,aa);

% Thus a first order filter is created with n = 2.
%[b,a] = butter(1, Fc/Fnyq);
%[zb,pb,kb] = butter(1, 2*pi*3e6, 's');
%[bb,ab] = zp2tf(zb,pb,kb);

% Plot the frequency response H(w):
%figure(1);
%w = linspace(0.4e6, 0.12e6)
%freqs(bb,ab, w);

% Remove DC component
%data = data - mean(data);
%data = filter(b,a,x);

%figure
%plot(data);

% number of data points
%N = length(data);

% numerical approx. of FT
%spec = fft(data);

% spacing between samples on freq. axis
%df = Fs/N

% min freq. for which fft is calculated
%min_f = -Fs/2;

% max freq. for which fft is calculated
%max_f = Fs/2 - df;

% horizontal values
%f = [min_f : df : max_f];

% should equal N
%size(f)

% magnitude of shifted spectrum
%y = abs(fftshift(spec));

% plot
%figure
%plot(f, y)
%xlabel('Frequency in [Hz]')
%ylabel('Amplitude')
%grid on

% Create a first-order Butterworth low pass.
% The returned vectors are of legth n.
% Thus a first order filter is created with n = 2.
%[b,a]=butter(1, Fc/(Fs/2));
%[zb,pb,kb] = butter(5, 2*pi*Fc, 's');
%[bb,ab] = zp2tf(zb,pb,kb);
%w = linspace(0.4e6, 0.6e6);
%figure(1);
%freqs(bb,ab,w);

% Compute its frequency response
%[H,w]=freqz(b, a, 512, Fsam);
%freqz_plot(w, abs(H));

% frequency response
%plot(w,abs(H));
%set(gca,'XScale','log')
%set(gca,'YScale','log')
%xlabel('Frequency (Hz)')
%ylabel('Filter Response')

% clear unused variables
%clear("Fnyq", "Fc");

% Clamp the sync tips to 0
%t=0:1/Fsam:5;
%input = repmat(x(1:end, 2), 5, 1); #-0.176;

% Apply the filter to the input signal and plot input and output.
%figure
%output=filter(b,a,data);
%plot(output);
%subplot(2, 1, 1)
%plot(filtered(:,1),";Impulse response;");
%subplot(2, 1, 2 )
%plot(transfer(:,1),";Butterworth Transfer;")
