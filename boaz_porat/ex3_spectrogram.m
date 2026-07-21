%try to reduce interpulations: tradeoffs and finding the most fitting
%window
%add time span and freq span

close all;

Fs = 2e3;
Fc = 50;
N = 1000;  
overlap = 0.25;

Ts = 1 / Fs;
dt = 1 / N;
RBW = Fs / N;

t = 0:Ts:Ts*(N-1);

signal = sin(2 * pi * Fc * t);

spec(signal, Fs, dt, RBW, overlap, hamming(RBW));


function spec(signal, Fs, dt, RBW, overlap, window)

N = Fs / RBW;

if N ~= 1 / dt
    error("RBW and dt must correspond")
end

Ts = 1 / Fs;
t = 0:Ts:Ts*(N-1) ;
f = -Fs/2:Fs/N:Fs/2 - Fs/N;
overlap_samples = floor(overlap * RBW);

buf_sig = buffer(signal, RBW, overlap_samples);
win_sig =  buf_sig .*  window;

fft_sig = fftshift(abs(fft(win_sig)));

un_buf_sig = fft_sig(:).';

figure;
plot(t,win_sig(:).');
figure;
plot(f, un_buf_sig);

figure;
imagesc(t, f, un_buf_sig);
figure;
spectrogram(win_sig(:).');
figure;
spectrogram(signal);

end


