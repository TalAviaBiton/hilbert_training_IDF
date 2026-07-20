fs = 10000;             % תדר דגימה (Hz)
duration_ms = 50;       % משך האות במילי-שניות
t = 0:1/fs:(duration_ms/1000 - 1/fs); 

f0 = 100;               % תדר התחלתי
f1 = 2000;              % תדר סופי

y = chirp(t, f0, duration_ms/1000, f1, 'linear');

bits = y > 0;

figure;
subplot(2,1,1);
plot(t*1000, y);
title('Chirp Signal (Waveform)');
xlabel('Time (ms)'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
spectrogram(y, 128, 120, 128, fs, 'yaxis');
title('Spectrogram of the Chirp Signal');