clc; clear; close all;

Fs = 2e3;
Fc = 50;
N = 1e3;  
overlap = 0.25;

Ts = 1 / Fs;
dt = 1 / N ;
RBW = N / 2;

t = 0:Ts:Ts*(N-1);

signal = sin(2 * pi * Fc * t);

t_start = 1; 
t_stop = N ; 
f_start = 1; 
f_stop = Fs / 2 ; 

spec(signal, Fs, dt, RBW, overlap, hamming(RBW), t_start, t_stop, f_start, f_stop);

% הנחתי שרוצים את האות כמו שהוא ולא בbb
function spec(signal, Fs, dt, RBW, overlap, window, t_start, t_stop, f_start, f_stop)

N = length(signal);
Ts = 1 / Fs;
t = 0:Ts:Ts * (N - 1) ;
f = -Fs / 2:Fs / N:Fs / 2 - Fs / N;
overlap_samples = floor(overlap * RBW);

if N ~= 1 / dt 
    error("RBW and dt must correspond to N");
end
if t_start < 1 || t_stop > N %could use some more validation tests but its matlab
    error("time span exeeded its bouds");
end
if f_start < 1 || f_stop > Fs / 2 %could use some more validation tests but its matlab
    error("frequncey span exeeded its bouds");
end

buf_sig = buffer(signal, RBW, overlap_samples);
win_sig =  buf_sig .*  window;

fft_sig = fftshift(abs(fft(win_sig, RBW, 1)));

un_buf_sig = fft_sig(:).';

ws = win_sig(:).';

wstp = ws(t_start:t_stop);
ubstp = un_buf_sig(f_start:f_stop);
stp = signal(t_start:t_stop);
ttp = t(t_start:t_stop);
ftp = f(f_start:f_stop);

figure;
plot(ttp,wstp);
figure;
plot(ftp, ubstp);

figure;
imagesc(ttp, ftp, ubstp);
figure;
spectrogram(wstp);
figure;
spectrogram(stp);

end

