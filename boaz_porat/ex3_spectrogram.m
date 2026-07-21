close all;

Fs = 1e3;
Fc = 50;
N = 1000;
overlap = 0;

Ts = 1 / Fs;
dt = 1 / N;
RBW = Fs / N;

t = 0:Ts:Ts*(N-1);

signal = sin(2 * pi * Fc * t);

spec(signal, Fs, dt, RBW, overlap, rectwin(RBW));


function spec(signal, Fs, dt, RBW, overlap, window)

N = Fs / RBW;

if N ~= 1 / dt
    error("RBW and dt must correspond")
end

Ts = 1 / Fs;
t = 0:Ts:Ts*(N-1) ;
f = -Fs/2:Fs/N:Fs/2 - Fs/N;

buf_sig = buffer(signal, RBW, overlap);
win_sig =  buf_sig .*  window;

fft_sig = fftshift(abs(fft(win_sig, RBW)));

figure
plot(t,(signal));
figure
plot(f, fft_sig);

figure
imagesc(t, f, fft_sig);
figure
spectrogram(signal);

end

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
