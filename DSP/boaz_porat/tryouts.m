close all;

fs = 50;             % קצב דגימה בהרץ
ts = 1/fs;
t = 0:ts:1;            % וקטור זמן למשך שנייה אחת
f = 5;                 % תדר אות הסינוס בהרץ
amplitude = 1;         % משרעת (אמפליטודה)

y = amplitude * sin(2 * pi * f * t); % חישוב האות

plot(t, y);            % הצגת הגרף
xlabel('זמן (שניות)');
ylabel('משרעת');
title('אות סינוס');
grid on;

N = length(y);
f = -fs/2:fs/N:fs/2 - fs/N;

fft_sig = fftshift(abs(fft(y)));

figure
plot(f, fft_sig);

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

%% 4

m = mindivisor(100, 3);

function m = mindivisor(N, a)
divs = N ./ (1:N);
for i = 1:N
    if mod(divs(i), 1) ~= 0
        divs(i) = inf(1,1);
    end
end
[~, idx] = min(abs(divs - a));
m= divs(idx);
end

%%

% הגדרות
L = 2;               % מקדם אינטרפולציה (Upsampling factor)
fs_new = 10000;      % תדר דגימה חדש לאחר האינטרפולציה
fc = 2500;           % תדר הקיטעון (התדר המקורי של האות, fs_old / 2)
N = 40;              % סדר המסנן

% נרמול תדר הקיטעון ביחס לתדר הדגימה החדש
Wn = fc / (fs_new / 2); % שווה ערך ל- 0.5

% תכנון המסנן (לעיתים מכפילים גם ב-L כדי לשמר את עוצמת האות)
b_interpolation = L * fir1(N, Wn, 'low');

% הצגת תגובת התדר
freqz(b_interpolation, 1, 1024, fs_new);
title('Interpolation Lowpass Filter');


%%
for i = 1:L
    left_extantion(i) = 2 * x(m1) - x(2 * m1 - i + m1);
end
% for i = m1 - L:m1
%     left_extantion(i) = 2 * x(m1) - x(2 * m1 - i);
% end

for i = 1:L
    right_extantion(i) = 2 * x(m2) - x(2 * m2 - i + m2);
end

