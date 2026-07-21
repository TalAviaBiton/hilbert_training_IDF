close all;

fs = 50;             % קצב דגימה בהרץ
ts = 1/fs
t = 0:ts:1;            % וקטור זמן למשך שנייה אחת
f = 5;                 % תדר אות הסינוס בהרץ
amplitude = 1;         % משרעת (אמפליטודה)

y = amplitude * sin(2 * pi * f * t); % חישוב האות

plot(t, y);            % הצגת הגרף
xlabel('זמן (שניות)');
ylabel('משרעת');
title('אות סינוס');
grid on;

function spec_try(signal, Fs, N)
Ts = 1 / Fs;
t = 0:Ts:Ts*(N-1) ;
f = -Fs/2:Fs/N:Fs/2 - Fs/N;

fft_sig = fftshift(abs(fft(signal)));

% figure
% plot(f, fft_sig);
% figure
% plot(t,(sig));
figure
imagesc(t, f, fft_sig);
figure
spectrogram(signal);
end