fc = 500;
Fs = 1e6; 
Wn = (2/Fs)*fc;

t = linspace(0,1,Fs);
x = cos(2*pi*1000*t) + cos(2*pi*150*t);

figure(1)
[b,a] = butter(1,Wn);
[h,f] = freqz(b,a,[],Fs);
plot(f,mag2db(abs(h)))
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid

figure(2)
y = filter(b,a,x);
plot(t,y,t,y)
xlim([0 0.1])
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original Signal','Filtered Data')