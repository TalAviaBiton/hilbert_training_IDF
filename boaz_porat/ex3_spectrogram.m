spec = create_spectrogram(signal, Fs, dt, RBW, overlap, symetryc_rect);

function spec = create_spectrogram(signal, Fs, dt, RBW, overlap, window)
% dt - required time resolution
% RBW - required frequency resolution
% overlap - percentage of overlap between two windows
% window - anonymous function that returns the window

f_start = 1e6;
f_stop = 2e6;
t_start = 0;
t_stop = 1e3;
Fc = 1.5e6;

%check RBW and dt correspond to each other


buf_sig= buffer(signal, RBW, overlap);

delta_f = Fs/RBW;
k = ceil(log2(delta_f));
N = 2 ^ k;
win = window(length(signal), N);

Ts = 1 / Fs;
t = t_start:1 / dt:t_stop - 1 / dt;
f = f_start:1 / RBW:f_stop - 1 / RBW;
imagesc(t, f, buf_sig);
set(buf_sig, 'ydir', 'normal')

end

function rect = symetryc_rect(len, num_of_samples)

pad = zeros(1,(len-num_of_samples)/2);
rect = [pad ones(1, num_of_samples) pad];

end