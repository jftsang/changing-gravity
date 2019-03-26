LyasFourier = fft(Lyas);
LyasFourier(1) = 0;

freqs = linspace(0, 2*pi/dt, nt);
periods = 2*pi ./ freqs;

figure;

subplot(2,1,1);
plot(ts, Lyas);

subplot(2,1,2);
plot(periods(1:floor(nt/2)), abs(LyasFourier)(1:floor(nt/2)));
xlim([0, 50]);
